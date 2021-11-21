# 1 better not minus, because centers should not be zero

# 2 newly added set_FirstRecBlock() cannot be adapted to fractures

# 3 newly added set_RegionMaxMinCorner().

In add fractures, only add_CircularFracture(), add_TriangularFracture, can adapt to it.

# 4 add fracture

不要太大，因为会生成一系列border points，根据region trim之后，不够三个border point就没有了

# 5 实现多个独立定义域分别切割的功能，且只留一个，或者用曙光搞这个功能也行

# 问题，切割完之后，添加完新顶点或者块体之后，表面polygon还是逆时针排列顶点么

这就很尴尬，切割之后的块体表面顶点排列不一定是逆时针。不过好在还连续，可是顺时针了就。搞了半天白搞，人家自己的搞得好好的。
