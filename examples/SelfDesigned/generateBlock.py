
from io import RawIOBase


class Block(object):
    def __init__(self, vertices, edges, faces, fixed=0) -> None:
        super().__init__()
        self.vertices = vertices  # coordinated of vertexes
        self.edges = edges  # v ids of edge nodes
        self.faces = faces  # v ids of face nodes
        self.fixed = fixed


def move(vertexes, offset):
    ni = len(vertexes)
    nj = len(offset)
    for i in range(ni):
        for j in range(nj):
            vertexes[i][j] += offset[j]
            # print(str(vertexes[i][j])+' ')
        # print('\n')
    return vertexes


def genBlk(Blocks: Block, s):
    nTotalBlocks = len(Blocks)
    nTotalVerts = 0
    nTotalEdges = 0
    nTotalFaces = 0
    nTotalPoints = 0
    for b in Blocks:
        nTotalVerts += len(b.vertices)
        nTotalEdges += len(b.edges)
        nTotalFaces += len(b.faces)
        for face in b.faces:
            nTotalPoints += len(face)
    f = open(s, 'w')
    f.write("#_DDA_DataFile_Version_1.0\n")
    f.write("POINTS_START_LENGTH "+str(nTotalBlocks)+'\n')
    auxId = 0
    for b in Blocks:
        f.write(str(auxId)+" "+str(len(b.vertices))+'\n')
        auxId += len(b.vertices)

    f.write(' '+'\n'+"POINTS "+str(nTotalVerts)+" double"+'\n')
    for b in Blocks:
        for v in b.vertices:
            f.write(str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')

    faces = []
    Blockfaces = []
    f.write(' '+'\n'+"FACES_NODELIST "+str(nTotalPoints)+'\n')
    auxId = 0
    auxId_f = 0
    auxId_blkf = 0
    for b in Blocks:
        Blockfaces.append([auxId_blkf, len(b.faces)])
        auxId_blkf += len(b.faces)

        for p in b.faces:
            faces.append([auxId_f, len(p)])
            auxId_f += len(p)
            for vId in p:
                f.write(str(vId)+' ')
        f.write('\n')

    f.write(' '+'\n'+"FACES "+str(nTotalFaces)+'\n')
    for face in faces:
        f.write(str(face[0])+' '+str(face[1])+' '+'\n')

    f.write(' '+'\n'+"BLOCK_FACES "+str(nTotalBlocks)+'\n')
    for bf in Blockfaces:
        f.write(str(bf[0])+' '+str(bf[1])+' '+'\n')

    f.write(' '+'\n'+"EDGES "+str(nTotalEdges)+'\n')
    Blockedges = []
    auxId = 0
    for b in Blocks:
        nEdges = len(b.edges)
        Blockedges.append([auxId, nEdges])
        for e in b.edges:
            f.write(str(e[0])+' '+str(e[1])+' ')
        auxId += nEdges
        f.write('\n')

    f.write("EDGES "+str(nTotalBlocks)+'\n')
    for e in Blockedges:
        f.write(str(e[0])+' '+str(e[1])+" \n")
    f.write(" \n")

    f.write("SCALARS aspectRatio float\n")
    for b in Blocks:
        f.write(str(3)+'\n')
    f.write(' '+'\n'+"SCALARS inscribedSphereRadius float"+'\n')

    for b in Blocks:
        f.write(str(1)+'\n')

    f.write(" \n"+"SCALARS fixedFlag int"+'\n')
    for b in Blocks:
        f.write(str(b.fixed)+'\n')

    f.write("SCALARS volume float\n")
    for b in Blocks:
        f.write("3\n")

    f.close()
