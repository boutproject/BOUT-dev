This script was modified from B. Dudson's original gridue_to_bout.py script. It creates a BOUT grid file from a GRIDUE file. 

It can be ran from the terminal.

To run it directly after the creation of a gridue file on INGRID the following code needs to be added to ingrid.py: 

```python
    def ExportBOUTgrid(self, gridue_file, bout_grid_name: str = 'bout_from_in.grd.nc', plotting = True, verbose = True, ignore_checks = False):
        """
        Export a BOUT grid file for the created grid.

        Parameters
        ----------
        gridue_file : str, mandatory
            Name of gridue file to convert to BOUT grid.

        bout_grid_name : str, optional
            Name of BOUT grid file to save.

        plotting : bool, optional
            If True, plot the gridue file before conversion.

        verbose : bool, optional
            If True, print verbose output during conversion.

        ignore_checks : bool, optional
            If True, ignore checks for gridue file format and structure.

        """
        bout_grid_name = gridue_file + "_" + bout_grid_name
        Convert_grids(gridue_file,bout_grid_name,plotting, verbose, ignore_checks)

```

And modifying the `ExportGridue` function to: 

```python
    if type(self.CurrentTopology) in [SNL]:
            if self.WriteGridueSNL(self.CurrentTopology.gridue_settings, fname):
                print(f"# Successfully saved gridue file as '{fname}'")
                self.ExportBOUTgrid(fname, 'bout_from_in.grd.nc', plotting = True, verbose = True, ignore_checks = False)
    elif type(self.CurrentTopology) in [SF15, SF45, SF75, SF105, SF135, SF165, UDN, CDN]:
        if self.WriteGridueDNL(self.CurrentTopology.gridue_settings, fname):
            print(f"# Successfully saved gridue file as '{fname}'")
            self.ExportBOUTgrid(fname, 'bout_from_in.grd.nc', plotting = True, verbose = True, ignore_checks = False)
```

*Note:* The created file will have 
`dimensions:
        x ;
        y ;
        z ;
        x2 ;
        y2 ;
        t = UNLIMITED ; //`

Where x has guard cells following INGRID's convention, and y doesn't have any. The converter takes the BOUT coordinates to be x -> y2 (and removes the guard cells) and y -> x2 (and adds the guard cells). So BOUT ends up using x2 as X and y2 as Y coordinates.