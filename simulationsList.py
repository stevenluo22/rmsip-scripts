def simulationsList():
    simulationslist = {}
    for r in range(1,17):
        simulationslist[f"run {r}"] = {
            "pdb": f"run_{r}/native.pdb", 
            "dcd": f"run_{r}/movie.dcd"
            }
    return simulationslist