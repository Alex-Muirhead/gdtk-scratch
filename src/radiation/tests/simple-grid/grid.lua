-- a simple grid

a = Vector3:new{x=0.0, y=0.0};
b = Vector3:new{x=1.0, y=0.0};
c = Vector3:new{x=1.0, y=10.0};
d = Vector3:new{x=0.0, y=10.0};

resolution = 27;

patch = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d};
grid = StructuredGrid:new{psurface=patch, niv=resolution+1, njv=resolution*10};

registerFluidGrid{grid=grid, fsTag="initial"};
