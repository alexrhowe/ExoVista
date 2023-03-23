/*  Bulirsch-Stoer integrator  

    (Unfortunately this is nearly a direct translation of a 
     FORTRAN code.  The original FORTRAN was ably written by
     Michael Henon.
     Translated to C by Jack Wisdom, adapted to C++/Cython by Alex Howe)

    must define EQUATIONS_OF_MOTION, MAX_NUMBER_EQUATIONS,
    INTEGRATION_EPSILON, number_equations,
    and possibly NAME in calling routine

*/

#include "specincld.h"
#include "Integrator.h"

Integrator::Integrator(vector<double> GMin, vector<double> stin, double tin, double dtin)
{
  GM = GMin;
  state = stin;
  t = tin;
  official_delta_t = dtin;
  
  n_equations = state.size();
  n_planets = n_equations / 6;
  epsilon = 1.e-12;
  max_n_planets = 30;
  max_n_equations = 6*max_n_planets;
  deriv = vector<double>(n_equations,0);
}


void Integrator::integrate()
{
  vector<double> local_state = state; //creates local versions of state to manipulate
  vector<double> state_deriv;
  static int lt_BS[10] = {1,2,3,4,6,8,12,16,24,32};
  vector<vector<double> > temp_BS(max_n_equations, vector<double>(12));
  vector<double> d_BS(6);
  double t_BS, delta_t_BS;

  int i_BS, m_BS, m1_BS, k_BS, i1max_BS, i1_BS;
  double xa_BS, xb_BS, varm_BS, h_BS, hd_BS, eta2_BS, dta_BS;
  double yb_BS, c_BS, b1_BS, den_BS, dtn_BS, b_BS, var_BS, varma_BS;
  double foo_BS, flt_BS;
  bool done = false;

  /* reset state variables--local_state set equal to state */
  t_BS = t;
  delta_t_BS = official_delta_t;

  xa_BS = t_BS;

  state_deriv = equations_of_motion(local_state);
 
  for(i_BS=0; i_BS<n_equations; i_BS++)
    {
      temp_BS[i_BS][1] = abs(local_state[i_BS]);
      if(temp_BS[i_BS][1] < epsilon)
	temp_BS[i_BS][1] = epsilon;
      temp_BS[i_BS][4] = state_deriv[i_BS];
      temp_BS[i_BS][0] = local_state[i_BS];
    }

  while(done==false)
    {
      xb_BS = delta_t_BS + xa_BS;
      for(m_BS=0; m_BS<10; m_BS++)
	{
	  flt_BS = lt_BS[m_BS];
	  varm_BS = 0.0;
	  m1_BS = min(m_BS, 6);
	  if(m1_BS != 0)
	    for(k_BS=0; k_BS<m1_BS; k_BS++)
	      {
		d_BS[k_BS] = flt_BS/lt_BS[m_BS-k_BS-1];
		d_BS[k_BS] *= d_BS[k_BS];
	      }
	  
	  h_BS = delta_t_BS / flt_BS;
	  hd_BS = 0.5 * h_BS;
	  for(i_BS=0; i_BS<n_equations; i_BS++)
	    {
	      temp_BS[i_BS][3] = temp_BS[i_BS][0];
	      local_state[i_BS] = temp_BS[i_BS][0] + hd_BS*temp_BS[i_BS][4];
	    }
	  
	  i1max_BS = 2*flt_BS - 1;
	  t_BS = xa_BS;
	  for(i1_BS=0; i1_BS<i1max_BS; i1_BS++)
	    {
	      t_BS += hd_BS;

	      state_deriv = equations_of_motion(local_state);
	      
	      for(i_BS=0; i_BS<n_equations; i_BS++)
		{
		  foo_BS = local_state[i_BS];
		  temp_BS[i_BS][1] = max(temp_BS[i_BS][1], abs(foo_BS));
		  if((temp_BS[i_BS][1]) == (0.0))
		    {
		      printf("FAIL: temp_BS[i_BS][1] = 0\n");
		      return;
		    }
		  eta2_BS = temp_BS[i_BS][3] + h_BS*state_deriv[i_BS];
		  temp_BS[i_BS][3] = local_state[i_BS];
		  local_state[i_BS] = eta2_BS;
		}
	    }

	  t_BS = xb_BS;
	  state_deriv = equations_of_motion(local_state);
	  
	  for(i_BS=0; i_BS<n_equations; i_BS++)
	    {
	      dta_BS = temp_BS[i_BS][11];
	      yb_BS = (temp_BS[i_BS][3] + local_state[i_BS] 
		       + hd_BS*state_deriv[i_BS]) / 2.0;
	      c_BS = yb_BS;
	      temp_BS[i_BS][11] = yb_BS;
	      if(m1_BS != 0)
		{
		  for(k_BS=0; k_BS<m1_BS; k_BS++)
		    {
		      b1_BS = d_BS[k_BS] * dta_BS;
		      den_BS = b1_BS - c_BS;
		      dtn_BS = dta_BS;
		      if(den_BS != 0)
			{
			  b_BS = (c_BS - dta_BS) / den_BS;
			  dtn_BS = c_BS * b_BS;
			  c_BS = b1_BS * b_BS;
			}
		      dta_BS = temp_BS[i_BS][11-k_BS-1];
		      temp_BS[i_BS][11-k_BS-1] = dtn_BS;
		      yb_BS += dtn_BS;
		    }
		  var_BS = abs(temp_BS[i_BS][2] - yb_BS) / temp_BS[i_BS][1];
		  varm_BS = max(varm_BS, var_BS);
		}
	      temp_BS[i_BS][2] = yb_BS;
	    }
	  if(m_BS < 3)
	    varma_BS = varm_BS; 
	  else if(varm_BS <= epsilon){ /* check for convergence or divergence */
	    done = true;
	    break;
	  }
	  else if(varm_BS >= varma_BS)
	    break;
	}
      if(done==false) delta_t_BS /= 2.0;
    }
  
  t_BS = xb_BS;
  for(i_BS=0; i_BS<n_equations; i_BS++)
    local_state[i_BS] = temp_BS[i_BS][2];
  
  /* reset state variables */
  t = t_BS;
  state = local_state;
  official_delta_t = delta_t_BS * 1.5 * pow(0.6, (double)(m_BS-m1_BS));
  
  return;
}


vector<double> Integrator::getState(){
  return state;
}

double Integrator::getTime(){
  return t;
}

double Integrator::getDeltaT(){
  return official_delta_t;
}

void Integrator::setDeltaT(double newdeltat){
  official_delta_t = newdeltat;
  return ;
}


vector<double> Integrator::equations_of_motion(vector<double> local_state)
{
  int c1, c2;
  double dx, dy, dz, r2, fij;
  vector<double> local_deriv(n_equations);

  // velocity from state
  for(int i=0; i<n_planets; i++){
    c1 = 6*i;
    local_deriv[c1]   = local_state[c1+3];
    local_deriv[c1+1] = local_state[c1+4];
    local_deriv[c1+2] = local_state[c1+5];
    local_deriv[c1+3] = 0.;
    local_deriv[c1+4] = 0.;
    local_deriv[c1+5] = 0.;
  }

  // pairwise accelerations
  for(int i=0; i<n_planets; i++){
    c1 = i*6;
    for(int j=0; j<i; j++){
      c2 = j*6;
      dx = local_state[c1]   - local_state[c2];
      dy = local_state[c1+1] - local_state[c2+1];
      dz = local_state[c1+2] - local_state[c2+2];
      r2 = dx*dx + dy*dy + dz*dz;
      fij = 1./pow(r2,1.5);
      local_deriv[c1+3] += -GM[j]*fij*dx;
      local_deriv[c1+4] += -GM[j]*fij*dy;
      local_deriv[c1+5] += -GM[j]*fij*dz;
      local_deriv[c2+3] +=  GM[i]*fij*dx;
      local_deriv[c2+4] +=  GM[i]*fij*dy;
      local_deriv[c2+5] +=  GM[i]*fij*dz;
    }
  }
  return local_deriv;
}
