import numpy
from matplotlib import pyplot
from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model
from prepare_figure import single_frame

def make_plot(disk1, filename):
    
    # close any open plots
    pyplot.close()
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=14)
    pyplot.xlim(-300, 300)
    pyplot.ylim(-300, 300)

    pyplot.scatter(disk1.x.value_in(units.kpc), disk1.y.value_in(units.kpc), 
                   alpha=1, s=5, lw=0)

    pyplot.savefig('plots/' + filename)

def make_galaxy(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk):
    converter=nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    galaxy1 = new_galactics_model(n_halo,
                                  converter,
                                  #do_scale = True,
                                  bulge_number_of_particles=n_bulge,
                                  disk_number_of_particles=n_disk)


    return galaxy1, converter

def simulate(galaxy1, converter, n_bulge, n_halo, t_end):
    
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    dynamics = Gadget2(converter, number_of_workers=8)
    dynamics.parameters.epsilon_squared = (100 | units.parsec)**2
    set1 = dynamics.particles.add_particles(galaxy1)
    dynamics.particles.move_to_center()
    
    disk1 = set1[:n_bulge] # This is just the thin and thick disks #:n_halo]
    bulge = set1[n_bulge:n_halo] # This is just the bulge
    disk_and_bulge = set1[:n_halo] # This is everything except dark matter halo
    
    # make a plot for each subset
    make_plot(disk1, "disk_0myr")
    make_plot(bulge, "bulge_0myr")
    make_plot(disk_and_bulge, "galaxy_0myr")
    
    sentinel = 1
    
    # use a while loop to simulate in parts, to get
    # snapshots of different parts of the simulation
    while sentinel*100 <= t_end.value_in(units.Myr):
        
        # evolve another 100 Myr
        dynamics.evolve_model(sentinel*100|units.Myr)
        
        # make and save the plots
        make_plot(disk1, "disk_" + str(100*sentinel) + "myr")
        make_plot(bulge, "bulge_" + str(100*sentinel) + "myr")
        make_plot(disk_and_bulge, "galaxy_" + str(100*sentinel) + "myr")
        
        # print to the terminal to see progress
        print('Done with t = ' + str(100*sentinel) + ' Myr')
        print()
        
        # add one to the sentinel
        sentinel += 1

    dynamics.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="M_galaxy", default = 1.0e12 | units.MSun,
                      help="Galaxy mass [%default]")
    result.add_option("-R", unit=units.kpc,
                      dest="R_galaxy", default = 10 | units.kpc,
                      help="Galaxy size [%default]")
    result.add_option("--n_bulge", dest="n_bulge", default = 10000,
                      help="number of stars in the bulge [%default]")
    result.add_option("--n_disk", dest="n_disk", default = 10000,
                      help="number of stars in the disk [%default]")
    result.add_option("--n_halo", dest="n_halo", default = 20000,
                      help="number of stars in the halo [%default]")
    result.add_option("--t_end", unit=units.Myr,
                      dest="t_end", default = 2000|units.Myr,
                      help="End of the simulation [%default]")
    return result

if __name__ == '__main__':
    o, arguments  = new_option_parser().parse_args()
    galaxy1, converter = make_galaxy(o.M_galaxy, o.R_galaxy,
                                                o.n_halo, o.n_bulge, o.n_disk)
    
    simulate(galaxy1, converter, int(o.n_bulge), int(o.n_halo), o.t_end)

