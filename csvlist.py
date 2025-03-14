def csvList():
    csvlist = {}
    for r in range(1,17):
        csvlist[f"run {r}"] = "rmsip.csv"
    return csvlist