from dataclasses import dataclass, field

@dataclass
class Settings():

    # imaging and spectroscopy parameters
    pixscale:     float = 0.002
    iwa:          float = 0.015
    iwa_tol:      float = 0.1
    npix:           int = 250
    specres:      float = 300.
    specrdisk:    float = 10.
    lambdamin:    float = 0.3
    lambdamax:    float = 1.0

    # planet population parameters
    seed:           int = None       # RNG seed
    emin:         float = 0.
    emax:         float = 0.
    imin:         float = 0.
    imax:         float = 5.
    sysimin:      float = 0.
    sysimax:      float = 180.
    sysPAmin:     float = 0.
    sysPAmax:     float = 360.    
    eecprob:      float = 1.0     # Probability that a planet in the EEC bounding box is an EEC.

    # disk model parameters
    ncomponents:     int = 2
    diskoff:        bool = False   # set to True to skip calculating the disk spectrum for fast testing
    minsize:       float = 0.1     # minimum disk particle size
    maxsize:       float = 1000.   # maximum disk particle size
    rdust_blowout: float = 0.5     # this could be updated for each star, but that's kinda BS anyway
    tsublimate :   float = 1500.
  
    # disk profile parameters
    density_ratio: float = 5.      # all components have the same density w/in this factor
    stability_factor: float = 3.
    rinner_mass_threshold: float = 100.
    dror_min:      float = 0.05     # dust belts must be at least this wide
    dror_max:      float = 0.3
    hor_min:       float = 0.03
    hor_max:       float = 0.2
    r_min:          list = field(default_factory=lambda: [0.5,5.,50.])  # min circumstellar distance of each disk component
    r_max:          list = field(default_factory=lambda: [5.,50.,500.]) # max circumstellar distance of each disk component

    # scene parameters
    timemax:      float = 1.e-10
    dt:           float = 10./365.25
    output_dir:     str = 'output'

    def __post_init__(self):
        self.pixscale_mas = self.pixscale * 1000.
