# Load a series of PDB files, each with its own XTC trajectory in VMD

# Define list of PDB and XTC file pairs and associated colour
set files {
    TEMP
}

# Loop through each pair and load them in VMD
foreach pair $files {
    set pdb_file [lindex $pair 0]  ;# Get PDB filename
    #set xtc_file [lindex $pair 1]  ;# Get XTC filename
    set pose_col [lindex $pair 1]  ;# Get VMD-Colour
    # Load the PDB file and get the molecule ID
    set molID [mol new $pdb_file type pdb waitfor all]
        mol delrep molID top
        mol selection "name BB and resid 1 to 439"
        mol representation Licorice 0.2 30
        mol color ColorID 2
        mol material AOChalky
        mol addrep top
        mol selection "name BB and resid 440 to 867"
        mol representation Licorice 0.2 30
        mol color ColorID 8
        mol material AOChalky
        mol addrep top
        mol selection "resname LIGNAME"
        mol representation Licorice 0.2 30
        mol color Name
        mol material AOShiny
        mol addrep top
        mol selection "(not resname LIGNAME) and within 7 of resname LIGNAME"
        mol representation Licorice 0.3 30
        mol color Name
        mol material AOShiny
        mol addrep top
        mol selection "resname LIGNAME"
        mol representation Licorice 0.3 30
        mol color ColorID $pose_col
        mol material AOShiny
        mol addrep top
    # # Load the corresponding XTC trajectory
    # mol addfile $xtc_file type xtc waitfor all molid $molID

    puts "Loaded $pdb_file"
}

puts "All structures and trajectories loaded successfully!"
