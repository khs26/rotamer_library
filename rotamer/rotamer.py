import sys

def rotamer_init(topology_filename):
    pass

if __name__ == '__main__':
    if sys.argv[1] == 'init':
        initialise(sys.argv[2])
    elif sys.argv[1] == 'step':
        take_step(sys.argv[2])