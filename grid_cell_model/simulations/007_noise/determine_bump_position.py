from models.place_input import PlaceCellInput

# WARNING: Check this values with the actual code in
# models.gc_net_nest.createGenericPlaceCells()
Ne_x = 34
Ne_y = 30
arenaSize = 180.0
gridSep = 60
connStdDev = gridSep / 2. / 6.
gridCenter = [0., 0.]

pc_input = PlaceCellInput(Ne_x, Ne_y, arenaSize, gridSep, [.0, .0],
        fieldSigma=connStdDev)


