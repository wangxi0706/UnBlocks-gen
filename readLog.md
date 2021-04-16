
# 1 better not minus, because centers should not be zero

# 2 newly added set_FirstRecBlock() cannot be adapted to fractures 

# 3 newly added set_RegionMaxMinCorner(). 
In add fractures, only add_CircularFracture(), add_TriangularFracture, can adapt to it.

# 4 add fracture
不要太大，因为会生成一系列border points，根据region trim之后，不够三个border point就没有了