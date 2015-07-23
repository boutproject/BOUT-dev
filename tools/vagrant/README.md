This is a basic vagrant config which you can use to compile and run BOUT++ on your local machine. For more info on vagrant see https://www.vagrantup.com/.

To use it, cd into this directory and run: 

	vagrant up                    # this creates the VM, installs BOUT including dependencies and runs all the tests. 
	vagrant ssh                   # logs in to the VM.

At this point you should be able to run one of the examples: 

	cd BOUT/examples/conduction   # find a simple example.
	make                          # compile the example.
	mpirun -np 2 ./conduction     # run the example.
