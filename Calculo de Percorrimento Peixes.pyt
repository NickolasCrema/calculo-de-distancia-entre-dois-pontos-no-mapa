# -*- coding: utf-8 -*-
#Version 1.0.0

import arcpy
import os
import copy
import sys


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = "toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        graphVerticesParam = arcpy.Parameter(
            name='graphVertices',
            displayName="Arquivo de vertices do grafo", 
            direction="Input", 
            parameterType="Required", 
            datatype="DEFeatureClass", 
        )

        graphEdgesParam = arcpy.Parameter(
            name='graphEdges',
            displayName='Arquivo de arestas do grafo',
            direction='Input',
            parameterType='Required',
            datatype='DEFeatureClass',
        )

        pointsFile = arcpy.Parameter(
            name='pointsFile',
            displayName='Arquivo de locais de recaptura',
            direction='Input',
            parameterType='Required',
            datatype='DEFeatureClass',
        )

        sqlExpression = arcpy.Parameter(
            name='sqlExpression',
            displayName='Filtro',
            direction='Input',
            parameterType='Optional',
            datatype='GPSQLExpression',
        )
        sqlExpression.parameterDependencies = [pointsFile.name]
        pointFile = arcpy.Parameter(
            name='pointFile',
            displayName='Arquivo de locais de soltura',
            direction='Input',
            parameterType='Required',
            datatype='DEFeatureClass',
        )

        checkBox = arcpy.Parameter(
            name='checkBox',
            displayName='Criar arquivo de trajetos',
            direction='Input',
            parameterType='Required',
            datatype='GPBoolean',
        )
        checkBox.value = False

        outputFileName = arcpy.Parameter(
            name='outputFileName',
            displayName='Nome destinado ao arquivo de trajetos',
            direction='Input',
            parameterType='Optional',
            datatype='GPString',
        )
        outputFileName.enabled = False

        outputPath = arcpy.Parameter(
            name='outputPath',
            displayName='Caminho destinado ao arquivo de trajetos',
            direction='Input',
            parameterType='Optional',
            datatype='DEWorkspace',
        )
        outputPath.enabled = False

        output = arcpy.Parameter(
            name='output',
            displayName='Arquivo de trajetos',
            direction='Input',
            parameterType='Optional',
            datatype='DEFeatureClass',
        )

        params = [graphVerticesParam, graphEdgesParam, pointsFile, sqlExpression, pointFile, checkBox, outputFileName, outputPath, output]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[5].value == True:
            parameters[8].enabled = False
            parameters[7].enabled = True
            parameters[6].enabled = True
        elif parameters[5].value == False:
            parameters[8].enabled = True
            parameters[7].enabled = False
            parameters[6].enabled = False
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        parameters = ParametersWrapper(parameters)
        graphVertices = parameters.graphVertices.valueAsText
        graphEdges = parameters.graphEdges.valueAsText
        srcPointFile = parameters.pointsFile.valueAsText
        destPointFile = parameters.pointFile.valueAsText
        outputFileName = parameters.outputFileName.valueAsText
        sqlExpression = parameters.sqlExpression.valueAsText
        output = parameters.output.valueAsText
        outputPath = parameters.outputPath.valueAsText
        checkBox = parameters.checkBox.value
        workSpace = os.path.dirname(graphEdges)
        fc = srcPointFile
        path = os.path.dirname(fc)
        #verifica se as features estão dentro de um dataset, caso estejam, volta um diretório
        if path[-4::] == '.gdb':
            editWorkspace = path
        else:
            editWorkspace = os.path.dirname('\\'.join(fc.split('\\')[:-1]))

        workspace = arcpy.env.scratchGDB

        spatialReference = arcpy.Describe(graphVertices).spatialReference
        inputSpatialReference = arcpy.Describe(srcPointFile).spatialReference

        #Verifica se o sistema de coordenadas dos inputs são iguais, caso contrário converte para o mesmo
        if inputSpatialReference != spatialReference:
            arcpy.management.Project(graphVertices, workspace + '\\vertices_project', inputSpatialReference)
            graphVertices = workspace + '\\vertices_project'
            
        spatialReference = arcpy.Describe(graphEdges).spatialReference
        if inputSpatialReference != spatialReference:
            arcpy.management.Project(graphEdges, workspace + '\\edges_project', inputSpatialReference)
            graphEdges = workspace + '\\edges_project'
            spatialReference = arcpy.Describe(graphEdges).spatialReference

        #Pega os nomes refentes ao OID@ e ao SHAPE@
        OIDSrc = arcpy.da.Describe(fc)['OIDFieldName']
        shapeSrc = arcpy.da.Describe(fc)['shapeFieldName']
        shapeDest = arcpy.da.Describe(destPointFile)['shapeFieldName']
        
        listVerticesFields = ['SHAPE@', 'Id', 'Edges', 'OID@']
        listEdgesFields = ['SHAPE@', 'Shape_Length', 'Id', 'Start', 'End', 'OID@']
        
        #dictionary para mapear o index do registro de vertice com base no seu OID
        dictGraphVerticesOID = {}
        #dictionary para mapear o index do registro de aresta com base no seu OID
        dictGraphEdgesOID = {}

        edgeCursor = arcpy.da.SearchCursor(graphEdges, listEdgesFields)
        vertexCursor = arcpy.da.SearchCursor(graphVertices, listVerticesFields)
        vertices = []
        edges = []

        for row in edgeCursor:
            edge = Graph.Edge(row[0], row[1], row[2], row[3], row[4])
            edges.append(edge)
            dictGraphEdgesOID[row[5]] = row[2]

        for row in vertexCursor:
            vertice = Graph.Vertex(row[0], row[1])
            vertice.fillEdges(row[2])
            vertices.append(vertice)
            dictGraphVerticesOID[row[3]] = row[1]

        graph = Graph()
        #popula a estrutura do grafo com seus dados
        graph.fillGraph(vertices, edges)
                
        #lista com os campos da tabela de recaptura
        listSrcFields = list(map(lambda field: field.name, arcpy.ListFields(srcPointFile)))
        #lista com os campos da tabela de soltura
        listDestFields = list(map(lambda field: field.name, arcpy.ListFields(destPointFile)))

        #dictionary para mapear o nome do campo com sua posição na tabela de recaptura
        dictSrc = {}
        i=0
        for field in listSrcFields:
            dictSrc[field] = i
            i+=1

        #dictionary para mapear o nome do campo com sua posição na tabela de soltura
        dictDest = {}
        i=0
        for field in listDestFields:
            dictDest[field] = i
            i+=1

        #registros do ponto de recaptura (com base na sqlExpression)
        srcCursor = arcpy.da.SearchCursor(srcPointFile, listSrcFields, where_clause=sqlExpression)
        #registros do ponto de soltura
        destCursor = arcpy.da.SearchCursor(destPointFile, listDestFields)

        #dictionary para mapear o shape do ponto de soltura com seu nome na tabela de recaptura
        dictLocaisSoltura = {}
        for row in destCursor:
            dictLocaisSoltura[row[dictDest['local_soltura']]] = row[dictDest[shapeDest]]

        # #verifica se o arquivo final de output já existe
        # if os.path.exists(outputPath + '\\' + outputFileName):
        #     output = outputPath + '\\' + outputFileName
        #     arcpy.AddMessage('existe')
        if checkBox == True:
            #cria o arquivo final de output
            arcpy.AddMessage('n existe')
            arcpy.AddMessage(outputPath + '\\' + outputFileName)
            output = arcpy.management.CreateFeatureclass(outputPath, outputFileName, 'POLYLINE', spatial_reference=spatialReference)
            arcpy.AddField_management(output, 'marca_externa_codigo', 'LONG', field_alias='Código da marca externa')
            arcpy.AddField_management(output, 'distancia_deslocamento_km', 'DOUBLE', field_alias='Distância de deslocamento do peixe (Km)')


        #percorre todos os registros (com base na sqlExpression) do arquivo de recaptura
        for row in srcCursor:
            if row[dictSrc['marca_externa_codigo']] == None:
                pass
            else:
                #ponto de recaptura
                source = row[dictSrc[shapeSrc]]
                #ponto de soltura
                dest = row[dictSrc['local_soltura']]
                arcpy.SetProgressor('default', "Carregando registro de OID: {0}...".format(row[dictSrc[OIDSrc]]))
                #verifica se o ponto de soltura é nulo, se for, é impossível calcular seu trajeto
                if dest != None:
                    destiny = dictLocaisSoltura[dest]
                    arcpy.SetProgressor('default', "Processando registro de OID: {0}...".format(row[dictSrc[OIDSrc]]))
                    #processa os dois pontos a serem calculados
                    res = ProcessPoints(graph, source, destiny, workspace, graphVertices, graphEdges, spatialReference, dictGraphEdgesOID, dictGraphVerticesOID)
                    #flag verifica o caso de undo/redo
                    flag = res[5]
                    #possíveis linhas geradas no Near do ponto do peixe com o grafo
                    polylinex = res[6]
                    polyliney = res[7]
                
                    #chama o algoritmo A* sobre os dois pontos
                    result = graph.A_Star(res[0], res[2])
                    
                    #cria feature auxiliar onde serão salvas as linhas do trajeto
                    outputEdges = arcpy.management.CreateFeatureclass(workSpace, 'outputedges', 'POLYLINE', spatial_reference=spatialReference)
                    #caso o resultado do algoritmo não seja nulo
                    if result != None:
                        arcpy.SetProgressor('default', "Escrevendo dados no registro de OID: {0}...".format(row[dictSrc[OIDSrc]]))
                        #adiciona as linhas percorridas na feature auxiliar
                        with arcpy.da.InsertCursor(outputEdges, ['SHAPE@', 'Shape_Length']) as cursorEdges:
                            
                            if polylinex != 0:
                                cursorEdges.insertRow([polylinex, polylinex.length])

                            for edge in result[1]:
                                cursorEdges.insertRow([graph.edges[edge].shape, graph.edges[edge].length])
                        
                            if polyliney != 0:
                                cursorEdges.insertRow([polyliney, polyliney.length])

                        #calcula a distância das linhas da feature de linhas percorridas em km
                        feature = arcpy.management.AddGeometryAttributes(outputEdges, 'LENGTH_GEODESIC', 'KILOMETERS')
                        codigo = row[dictSrc['marca_externa_codigo']]
                        resultCursor = arcpy.da.UpdateCursor(srcPointFile, listSrcFields, where_clause=sqlExpression)

                        #inicia uma sessão de edição no workspace
                        edit = arcpy.da.Editor(editWorkspace)
                        edit.startEditing(False, False)
                        edit.startOperation()
                        length = 0
                        lines = []
                        #corrige o arquivo de linhas percorridas e soma suas distâncias
                        with arcpy.da.UpdateCursor(feature, ['LENGTH_GEO', 'SHAPE@']) as cursor:
                            for rowAux in cursor:
                                if rowAux[1] != None:
                                    lines.append(rowAux[1])
                                if rowAux[0] == None:
                                    cursor.deleteRow()
                                else:
                                    length += rowAux[0]

                        unsplitedLines = joinShapes(lines)
                        #insere as linhas percorridas no arquivo final de output
                        with arcpy.da.InsertCursor(output, ['SHAPE@', 'marca_externa_codigo', 'distancia_deslocamento_km']) as cursor:
                            cursor.insertRow([unsplitedLines, codigo, length])
                        #modifica a distância percorrida no arquivo de soltura
                        for rowR in resultCursor:
                            if rowR[dictSrc['marca_externa_codigo']] == codigo:
                                rowR[dictSrc['distancia_deslocamento_km']] = length
                                resultCursor.updateRow(rowR)
                                break
                            else:
                                continue
                        
                        edit.stopOperation ()
                        edit.stopEditing (True)
                        
                    #volta o grafo para o estado em que estava antes de realizar o processamento
                        if flag == 0:
                            arcpy.SetProgressor('default', "Registro de OID: {0} processado com sucesso...".format(row[dictSrc[OIDSrc]]))
                            arcpy.AddMessage('Registro de OID: ' + str(row[dictSrc[OIDSrc]]) + ' efetuado com sucesso!')
                            continue
                        elif flag == 1:
                            undo = res[4][0]
                            start = undo[0]
                            ed = undo[1]
                            end = undo[2]
                            for edge in res[4][1]:
                                graph.edges.pop(edge)
                            for vert in res[4][2]:
                                graph.vertices.pop(vert)
                            graph.vertices[start.index].renewVertex(copy.deepcopy(start))
                            graph.edges[ed.index].renewEdge(copy.deepcopy(ed))
                            graph.vertices[end.index].renewVertex(copy.deepcopy(end))
                            del start
                            del ed
                            del end
                        elif flag == 2:
                            undo = res[4][0]
                            undo2 = res[4][3]
                            start1 = undo[0]
                            ed1 = undo[1]
                            end1 = undo[2]
                            start2 = undo2[0]
                            ed2 = undo2[1]
                            end2 = undo2[2]
                            for edge in res[4][1]:
                                graph.edges.pop(edge)
                            for edge in res[4][4]:
                                graph.edges.pop(edge)
                            for vert in res[4][2]:
                                graph.vertices.pop(vert)
                            for vert in res[4][5]:
                                graph.vertices.pop(vert)
                            graph.vertices[start1.index].renewVertex(copy.deepcopy(start1))
                            graph.edges[ed1.index].renewEdge(copy.deepcopy(ed1))
                            graph.vertices[end1.index].renewVertex(copy.deepcopy(end1))
                            graph.vertices[start2.index].renewVertex(copy.deepcopy(start2))
                            graph.edges[ed2.index].renewEdge(copy.deepcopy(ed2))
                            graph.vertices[end2.index].renewVertex(copy.deepcopy(end2))
                            del start1
                            del ed1
                            del end1
                            del start2
                            del ed2
                            del end2
                        elif flag == 4:
                            undo = res[4][0]
                            start = undo[0]
                            ed = undo[1]
                            end = undo[2]
                            for edge in res[4][1]:
                                graph.edges.pop(edge)
                            for vert in res[4][2]:
                                graph.vertices.pop(vert)
                            graph.vertices[start.index].renewVertex(copy.deepcopy(start))
                            graph.edges[ed.index].renewEdge(copy.deepcopy(ed))
                            graph.vertices[end.index].renewVertex(copy.deepcopy(end))
                            del start
                            del ed
                            del end
                        arcpy.SetProgressor('default', "Registro de OID: {0} processado com sucesso...".format(row[dictSrc[OIDSrc]]))
                        arcpy.AddMessage('Registro de OID: ' + str(row[dictSrc[OIDSrc]]) + ' efetuado com sucesso!')
                        del edit
                        arcpy.management.Delete(outputEdges)
                        arcpy.management.Delete(unsplitedLines)
                        arcpy.management.Delete(feature)
                    else:
                        arcpy.AddMessage('Registro de OID: ' + str(row[dictSrc[OIDSrc]]) + ' falhou!')
                        continue
                
        arcpy.management.Delete(graphVertices)
        arcpy.management.Delete(graphEdges)
                
        return



def pointSqrDistance(A, B):
    return (A.X - B.X)**2 + (A.Y - B.Y)**2

def heuristic(a, b):
    return pointSqrDistance(a.point(), b)**0.5

#constrói o caminho percorrido no algoritmo A*
#Parametros:
#   passedThrough: set com o id das arestas percorridas
#   cameFrom: set com o id dos vertices percorridos
#   current: vertice atual
#   start: vertice inicial
#retorna uma tupla;
#   total_path_vertex: lista com o id dos vertices percorridos
#   total_path_edge: lista com o id das arestas percorridas
def reconstruct_path(passedThrough, cameFrom, current, start):
    total_path_vertex = [current]
    total_path_edge = []
    while current in cameFrom:
        aux = passedThrough[current]
        current = cameFrom[current]
        total_path_vertex.insert(0, current)
        total_path_edge.insert(0, aux)
        if current == start:
            break
        
        
    return total_path_vertex, total_path_edge

def distanceSqr(A, B):
    return (A.X - B.X)**2 + (A.Y - B.Y)**2

def joinParts(shape):
    if shape.partCount == 1:
        return shape.getPart(0)
    #Precisa clonar?
    joinedParts = shape.getPart(0)
    for i in range(1, shape.partCount):
        joinedParts.extend(shape.getPart(i))
    return joinedParts

def joinShapes(shapes):
    if len(shapes) == 1:
        return shapes[0]
    #Discover proper orientation
    parts0 = joinParts(shapes[0])
    parts1 = joinParts(shapes[0])
    A0 = parts0[0]
    B0 = parts0[1]
    A1 = parts1[0]
    B1 = parts1[1]
    dA = min(distanceSqr(A0, A1), distanceSqr(A0, B1))
    dB = min(distanceSqr(B0, A1), distanceSqr(B0, B1))
    if dA < dB:
        allParts = arcpy.Array(reversed(parts0))
    else:
        allParts = parts0
    for shape in shapes[1:]:
        parts1 = joinParts(shape)
        A1 = parts1[0]
        B1 = parts1[1]
        B0 = allParts[-1]
        if distanceSqr(B0, A1) <= distanceSqr(B0, B1):
            allParts.extend(parts1)
        else:
            allParts.extend(reversed(parts1))
    spatialReference = shapes[0].spatialReference
    return arcpy.Polyline(allParts, spatial_reference = spatialReference)

#estrutura de grafo
class Graph:
    #estrutura de aresta
    class Edge:
        def __init__(self, shape, length, index, start, end):
            self.shape = shape
            self.length = length
            self.start = start
            self.end = end
            self.index = index
        
        def __repr__(self):
            return f"({self.index}, {self.start} -> {self.end})"
            
        def renewEdge(self, edge):
            self.shape = edge.shape
            self.length = edge.length
            self.start = edge.start
            self.end = edge.end
    #estrutura de vertice
    class Vertex:
        def __init__(self, shape, index):
            self.shape = shape
            self.index = index
            self.edges = []
        
        def fillEdges(self, edges):
            edgeSplit = edges.split(';')
            for edge in edgeSplit:
                self.edges.append(edge)        
        
        def point(self):
            return self.shape.getPart(0)
        
        def __repr__(self):
            aux = ", ".join(str(i) for i in self.edges)
            return f"({self.index}, [{aux}])"
        
        def renewVertex(self, vertex):
            self.shape = vertex.shape
            self.edges = vertex.edges

    def __init__ (self):
        self.vertices = {}
        self.edges = {}

    def fillGraph(self, vertices, edges):
        for vertex in vertices:
            self.vertices[vertex.index] = vertex
        for edge in edges:
            self.edges[edge.index] = edge
        
        
    def __repr__(self):
        return f"{{{self.Edge}, {self.Vertex}}}"

    #A* encontra o menor caminho entre o ponto inicial e final
    #Parametros:
    #   start: estrutura de vertice do ponto inicial
    #   goal: estrutura de vertice do ponto final
    #retorno void
    def A_Star(self, start, goal):
        goalPoint = self.vertices[goal.index].point()
        def h(nodePoint):
            return heuristic(nodePoint, goalPoint)
        # The set of discovered nodes that may need to be (re-)expanded.
        # Initially, only the start node is known.
        # This is usually implemented as a min-heap or priority queue rather than a hash-set.
        openSet = {start.index}

        alreadyPassed = {}
        for key in self.edges:
            alreadyPassed[key] = False
        
        # For node n, cameFrom[n] is the node immediately preceding it on the cheapest path from start
        # to n currently known.
        cameFrom = {}

        passedThrough = {}
        # cameFrom[start]

        # For node n, gScore[n] is the cost of the cheapest path from start to n currently known.
        gScore = {} #if node not on map, the default value is infinity
                    #so, use as gScore.get(node, float('inf'))
        gScore[start.index] = 0

        # For node n, fScore[n] := gScore[n] + h(n). fScore[n] represents our current best guess as to
        # how short a path from start to finish can be if it goes through n.
        fScore = {} #map with default value of Infinity
        fScore[start.index] = h(start)

        while openSet:
            # This operation can occur in O(Log(N)) time if openSet is a min-heap or a priority queue
            current = min(openSet, key=lambda node: fScore[node]) #the node in openSet having the lowest fScore[] value
    
            aux = self.vertices[current]
            if goal == self.vertices[current]:
                return reconstruct_path(passedThrough, cameFrom, current, start)

            openSet.remove(current)

            current = self.vertices[current]
            for neighborPath in current.edges:
                # d(current,neighbor) is the weight of the edge from current to neighbor
                # tentative_gScore is the distance from start to the neighbor through current
                tentative_gScore = gScore[current.index] + self.edges[int(neighborPath)].length

                nextVert = 0
                if self.edges[int(neighborPath)].start == current.index:
                    nextVert = self.edges[int(neighborPath)].end
                elif self.edges[int(neighborPath)].end == current.index:
                    nextVert = self.edges[int(neighborPath)].start


                if nextVert not in gScore.keys():
                    gScore[int(nextVert)] = sys.maxsize

                if tentative_gScore < gScore[int(nextVert)]:
                    # This path to neighbor is better than any previous one. Record it!

                    gScore[int(nextVert)] = tentative_gScore
                    fScore[int(nextVert)] = tentative_gScore + h(self.vertices[nextVert])
                    if neighborPath not in openSet and alreadyPassed[int(neighborPath)] == False:
                       
                        passedThrough[int(nextVert)] = int(neighborPath)
                        cameFrom[int(nextVert)] = current.index
                        alreadyPassed[int(neighborPath)] = True
                        openSet.add(nextVert)
                        
        # Open set is empty but goal was never reached
        return
    

class ParametersWrapper(object):
    """Empacota os parâmetros para permitir acesso por nome."""
    def __init__(self, parameters):
        for p in parameters:
            self.__dict__[p.name] = p

#Aproxima os pontos x e y da malha e processa os dados de acordo com cada caso
#Parametros:
#   graph: estrutura do grafo
#   pointx: ponto x a ser aproximado da malha
#   pointy: ponto y a ser aproximado da malha
#   workspace: caminho do diretório fonte para salvar features auxiliares
#   verticesPath: caminho do arquivo de vértices
#   edgesPath: caminho do arquivo de arestas
#   spatialReference: sistema de coordenadas
#   dictE: dictionary que mapeia as arestas com base no OID
#   dictV: dictionary que mapeia os vertices com base no OID
#retorna uma tupla;
#   idx: o vértice referente ao ponto x
#   closestDistancex: a distancia do ponto x até o novo vértice da malha
#   idy: o vértice referente ao ponto y
#   closestDistancey: a distancia do ponto y até o novo vértice da malha
#   undo: uma lista de listas composta de pontos, arestas, id para reversão das alterações feitas no grafo,
#         na forma de [[pontoinicial, aresta, pontofinal], [id], [id]]
#   flag: uma variável de inteiro para definir qual caso de processamento foi necessário
def ProcessPoints(graph, pointx, pointy, workspace, verticesPath, edgesPath, spatialReference, dictE, dictV):
    pointFeaturex = arcpy.management.CreateFeatureclass(workspace, 'pointFeature', 'POINT', spatial_reference=spatialReference)
    pointFeaturey= arcpy.management.CreateFeatureclass(workspace, 'pointFeature2', 'POINT', spatial_reference=spatialReference)
    malha = [verticesPath, edgesPath]
    verticesName = verticesPath.split('\\')[-1]
    edgesName = edgesPath.split('\\')[-1]

    arrayx = arcpy.Array()
    arrayy = arcpy.Array()
    point = arcpy.Point()
    point.X = pointx[0]
    point.Y = pointx[1]
    arrayx.add(point)
    point.X = pointy[0]
    point.Y = pointy[1]
    arrayy.add(point)

    #joga o pontox em uma feature para processamento
    with arcpy.da.InsertCursor(pointFeaturex, ['SHAPE@']) as cursor:
        cursor.insertRow([pointx])
    
    #joga o pontoy em uma feature para processamento
    with arcpy.da.InsertCursor(pointFeaturey, ['SHAPE@']) as cursor:
        cursor.insertRow([pointy])
    
    #realiza o Near do pontox com a malha do grafo
    arcpy.analysis.Near(pointFeaturex, malha, location=True, method='GEODESIC')
    #realiza o Near do pontoy com a malha do grafo
    arcpy.analysis.Near(pointFeaturey, malha, location=True, method='GEODESIC')

    #pega as informações resultantes do pontox
    with arcpy.da.SearchCursor(pointFeaturex, ['NEAR_FID', 'NEAR_DIST', 'NEAR_X', 'NEAR_Y', 'NEAR_FC']) as cursor:
        for row in cursor:
            pointInfox = row

    #pega as informações resultantes do pontoy
    with arcpy.da.SearchCursor(pointFeaturey, ['NEAR_FID', 'NEAR_DIST', 'NEAR_X', 'NEAR_Y', 'NEAR_FC']) as cursor:
        for row in cursor:
            pointInfoy = row
    
    closestIDx = pointInfox[0]
    closestDistancex = pointInfox[1] * 0.0001
    closestXx = pointInfox[2]
    closestYx = pointInfox[3]
    closestTypex = pointInfox[4].split('\\')[-1]

    closestIDy = pointInfoy[0]
    closestDistancey = pointInfoy[1] * 0.0001
    closestXy = pointInfoy[2]
    closestYy = pointInfoy[3]
    closestTypey = pointInfoy[4].split('\\')[-1]

    pointx = arcpy.Point(closestXx, closestYx)
    pointx = arcpy.PointGeometry(pointx, spatial_reference= spatialReference)

    pointy = arcpy.Point(closestXy, closestYy)
    pointy = arcpy.PointGeometry(pointy, spatial_reference= spatialReference)

    polylinex = 0
    polyliney = 0
    idx = 0 
    idy = 0
    index = 0, 0, 0
    flag = 0
    undo = []

    if closestDistancex != 0:
        # arcpy.AddMessage('closestDistancex: ' + str(closestDistancex))
        point.X = closestXx
        point.Y = closestYx
        arrayx.add(point)
        polylinex = arcpy.Polyline(arrayx, spatial_reference= spatialReference)
    if closestDistancey != 0:
        # arcpy.AddMessage('closestDistancey: ' + str(closestDistancey))
        point.X = closestXy
        point.Y = closestYy
        arrayy.add(point)
        polyliney = arcpy.Polyline(arrayy, spatial_reference= spatialReference)

    #salva os dados da aresta x
    if closestTypex == edgesName:
        undox = [copy.deepcopy(graph.vertices[graph.edges[dictE[closestIDx]].start]), copy.deepcopy(graph.edges[dictE[closestIDx]]), copy.deepcopy(graph.vertices[graph.edges[dictE[closestIDx]].end])]
    
    #salva os dados da aresta y
    if closestTypey == edgesName:
        undoy = [copy.deepcopy(graph.vertices[graph.edges[dictE[closestIDy]].start]), copy.deepcopy(graph.edges[dictE[closestIDy]]), copy.deepcopy(graph.vertices[graph.edges[dictE[closestIDy]].end])]

    
    #primeiro caso, se os dois pontos coincidem na mesma aresta
    if closestIDx == closestIDy and closestTypex == closestTypey and closestTypex == edgesName:
        flag = 4
        index = TwoTimesSplitEdge(graph, dictE[closestIDx], pointx, pointy, workspace, spatialReference)
        idx = index[0]
        idy = index[1]
        undo.append(undox)
        undo.append(index[2])
        undo.append(index[3])
    
    else:
        #o ponto x coincide em uma aresta
        if closestTypex == edgesName:
            flag+= 1
            index = SplitEdge(graph, dictE[closestIDx], pointx, workspace, spatialReference)
            idx = index[0]
            undo.append(undox)
            undo.append(index[1])
            undo.append(index[2])

            #segundo caso, ambos os pontos coincidem em uma aresta, porém não a mesma
            if closestTypey == edgesName:
                flag+= 1
                index = SplitEdge(graph, dictE[closestIDy], pointy, workspace, spatialReference)
                idy = index[0]
                undo.append(undoy)
                undo.append(index[1])
                undo.append(index[2])
            
            #terceiro caso, só um ponto coincide em uma aresta, o outro já é um vertice no grafo
            elif closestTypey == verticesName:
                idy = graph.vertices[dictV[closestIDy]]
                flag+= 0

        #o ponto x já é um vertice no grafo
        elif closestTypex == verticesName:
            idx = graph.vertices[dictV[closestIDx]]
            flag+= 0

            #terceiro caso, só um ponto coincide em uma aresta, o outro já é um vertice no grafo
            if closestTypey == edgesName:
                flag+= 1
                index = SplitEdge(graph, dictE[closestIDy], pointy, workspace, spatialReference)
                idy = index[0]
                undo.append(undoy)
                undo.append(index[1])
                undo.append(index[2])
        
            #quarto caso, ambos os pontos não coincidem em uma aresta, ambos já são vertices no grafo.
            elif closestTypey == verticesName:
                idy = graph.vertices[dictV[closestIDy]]
                flag+= 0

    return idx, closestDistancex, idy, closestDistancey, undo, flag, polylinex, polyliney

#processa os dois pontos quando incidem na mesma aresta
#cria dois vértices em cima da mesma aresta e divide a aresta em três
#Parametros:
#   graph: estrutura do grafo
#   closestIDx: inteiro representando o id da aresta x mais próxima do grafo na malha obtida pela função Near
#   closestIDy: inteiro representando o id da aresta y mais próxima do grafo na malha obtida pela função Near
#   pointx: PointGeometry do ponto x
#   pointy: PointGeometry do ponto y
#   workspace: caminho do diretório fonte para salvar features auxiliares
#   spatialReference: sistema de coordenadas
#retorna uma tupla;
#   graph.vertices[vertx.index]: o vertice com o novo ponto x criado sobre a aresta
#   graph.vertices[verty.index]: o vertice com o novo ponto y criado sobre a aresta
#   undoNewEdge: uma lista com o id das novas arestas criadas, para reversão posterior
#   undoNewVert: uma lista com o id dos novos vertices criados, para reversão posterior
def TwoTimesSplitEdge(graph, closestIDx, pointx, pointy, workspace, spatialReference):
    edgex = graph.edges[closestIDx]

    Edgex = arcpy.management.CreateFeatureclass(workspace, "Edgex", 'POLYLINE', spatial_reference=spatialReference)
    shape = ['SHAPE@']
    with arcpy.da.InsertCursor(Edgex, shape) as cursor:
        cursor.insertRow([edgex.shape])
    singlePointx = arcpy.management.CreateFeatureclass(workspace, "singlePointx", 'POINT', spatial_reference=spatialReference)
    with arcpy.da.InsertCursor(singlePointx, shape) as cursor:
        cursor.insertRow([pointx])
    saidax = arcpy.management.CreateFeatureclass(workspace, 'saidax', 'POLYLINE', spatial_reference=spatialReference)

    arcpy.management.SplitLineAtPoint(Edgex, singlePointx, saidax)
    edgesx = []
    with arcpy.da.SearchCursor(saidax, ['SHAPE@', 'Shape_Length']) as cursor:
        for row in cursor:
            edgesx.append(row)

    vertx = Graph.Vertex(pointx, max([i for i in graph.vertices.keys()]) + 1)
    lastIndexEdgex = max([i for i in graph.edges.keys()]) + 1

    graph.vertices[vertx.index] = vertx
    graph.vertices[vertx.index].fillEdges(str(closestIDx) + ';' + str(lastIndexEdgex))

    graph.edges[lastIndexEdgex] = Graph.Edge(edgesx[0][0], edgesx[0][1], lastIndexEdgex, graph.edges[closestIDx].start, vertx.index)
    if str(closestIDx) in graph.vertices[graph.edges[closestIDx].start].edges:
        graph.vertices[graph.edges[closestIDx].start].edges.remove(str(closestIDx))
    if str(lastIndexEdgex) not in graph.vertices[graph.edges[closestIDx].start].edges:
        graph.vertices[graph.edges[closestIDx].start].edges.append(str(lastIndexEdgex))
    graph.edges[closestIDx] = Graph.Edge(edgesx[1][0], edgesx[1][1], closestIDx, vertx.index, graph.edges[closestIDx].end)
    
    edg1 = []
    edg2 = []

    for e in graph.edges[closestIDx].shape.getPart()[0]:
        edg1.append(str(e))

    for e in graph.edges[lastIndexEdgex].shape.getPart()[0]:
        edg2.append(str(e))

    if str(pointy.getPart()) in edg1:
        edgey = graph.edges[closestIDx]
        Edgey = arcpy.management.CreateFeatureclass(workspace, 'Edgey', 'POLYLINE', spatial_reference=spatialReference)
        with arcpy.da.InsertCursor(Edgey, shape) as cursor:
            cursor.insertRow([edgey.shape])
        singlePointy = arcpy.management.CreateFeatureclass(workspace, 'singlePointy', 'POINT', spatial_reference=spatialReference)
        with arcpy.da.InsertCursor(singlePointy, shape) as cursor:
            cursor.insertRow([pointy])
        saiday = arcpy.management.CreateFeatureclass(workspace, 'saiday', 'POLYLINE', spatial_reference=spatialReference)

        arcpy.management.SplitLineAtPoint(Edgey, singlePointy, saiday)
        edgesy = []
        with arcpy.da.SearchCursor(saiday, ['SHAPE@', 'Shape_Length']) as cursor:
            for row in cursor:
                edgesy.append(row)
        
        verty = Graph.Vertex(pointy, max([i for i in graph.vertices.keys()]) + 1)
        lastIndexEdgey = max([i for i in graph.edges.keys()]) + 1

        graph.vertices[verty.index] = verty
        graph.vertices[verty.index].fillEdges(str(closestIDx) + ';' + str(lastIndexEdgey))

        graph.edges[lastIndexEdgey] = Graph.Edge(edgesy[0][0], edgesy[0][1], lastIndexEdgey, graph.edges[closestIDx].start, verty.index)
        if str(closestIDx) in graph.vertices[graph.edges[closestIDx].start].edges:
            graph.vertices[graph.edges[closestIDx].start].edges.remove(str(closestIDx))
        if str(lastIndexEdgey) not in graph.vertices[graph.edges[closestIDx].start].edges:
            graph.vertices[graph.edges[closestIDx].start].edges.append(str(lastIndexEdgey))
        graph.edges[closestIDx] = Graph.Edge(edgesy[1][0], edgesy[1][1], closestIDx, verty.index, graph.edges[closestIDx].end)
    
    elif str(pointy.getPart()) in edg2:
        edgey = graph.edges[lastIndexEdgex]
        Edgey = arcpy.management.CreateFeatureclass(workspace, 'Edgey', 'POLYLINE', spatial_reference=spatialReference)
        with arcpy.da.InsertCursor(Edgey, shape) as cursor:
            cursor.insertRow([edgey.shape])
        singlePointy = arcpy.management.CreateFeatureclass(workspace, 'singlePointy', 'POINT', spatial_reference=spatialReference)
        with arcpy.da.InsertCursor(singlePointy, shape) as cursor:
            cursor.insertRow([pointy])
        saiday = arcpy.management.CreateFeatureclass(workspace, 'saiday', 'POLYLINE', spatial_reference=spatialReference)
        
        arcpy.management.SplitLineAtPoint(Edgey, singlePointy, saiday)
        edgesy = []
        with arcpy.da.SearchCursor(saiday, ['SHAPE@', 'Shape_Length']) as cursor:
            for row in cursor:
                edgesy.append(row)

        verty = Graph.Vertex(pointy, max([i for i in graph.vertices.keys()]) + 1)
        lastIndexEdgey = max([i for i in graph.edges.keys()]) + 1

        graph.vertices[verty.index] = verty
        graph.vertices[verty.index].fillEdges(str(lastIndexEdgex) + ';' + str(lastIndexEdgey))

        graph.edges[lastIndexEdgey] = Graph.Edge(edgesy[0][0], edgesy[0][1], lastIndexEdgey, graph.edges[lastIndexEdgex].start, verty.index)
        if str(lastIndexEdgex) in graph.vertices[graph.edges[lastIndexEdgex].start].edges:
            graph.vertices[graph.edges[lastIndexEdgex].start].edges.remove(str(lastIndexEdgex))
        if str(lastIndexEdgey) not in graph.vertices[graph.edges[closestIDx].start].edges:
            graph.vertices[graph.edges[lastIndexEdgex].start].edges.append(str(lastIndexEdgey))
        graph.edges[lastIndexEdgex] = Graph.Edge(edgesy[1][0], edgesy[1][1], lastIndexEdgex, verty.index, graph.edges[lastIndexEdgex].end)
    else:
        arcpy.AddMessage('Um erro inesperado aconteceu...')
    undoNewEdge = [lastIndexEdgex, lastIndexEdgey]
    undoNewVert = [vertx.index, verty.index]
    arcpy.management.Delete(Edgex)
    arcpy.management.Delete(singlePointx)
    arcpy.management.Delete(Edgey)
    arcpy.management.Delete(singlePointy)

    return graph.vertices[vertx.index], graph.vertices[verty.index], undoNewEdge, undoNewVert

#processa um ponto quando incide em uma aresta
#cria um vértice em cima da aresta e divide a aresta em duas
#Parametros:
#   graph: estrutura do grafo
#   closestID: inteiro representando o id da aresta mais próxima do grafo na malha obtida pela função Near
#   point: PointGeometry do ponto mais próximo do grafo na malha
#   workspace: caminho do diretório fonte para salvar features auxiliares
#   spatialReference: sistema de coordenadas
#retorna uma tupla;
#   graph.vertices[vert.index]: o vertice com o novo ponto criado sobre a aresta
#   undoNewEdge: uma lista com o id da nova aresta, para reversão posterior
#   undoNewVert: uma lista com o id do novo vertice, para reversão posterior
def SplitEdge(graph, closestID, point, workspace, spatialReference):
    edge = graph.edges[closestID]
    Edge = arcpy.management.CreateFeatureclass(workspace, "Edge", 'POLYLINE', spatial_reference=spatialReference)
    shape = ['SHAPE@']
    with arcpy.da.InsertCursor(Edge, shape) as cursor:
        cursor.insertRow([edge.shape])
    pointx = point
    singlePoint = arcpy.management.CreateFeatureclass(workspace, 'singlePoint', 'POINT', spatial_reference=spatialReference)
    with arcpy.da.InsertCursor(singlePoint, shape) as cursor:
        cursor.insertRow([pointx])

    saida = str(arcpy.management.CreateFeatureclass(workspace, 'saida', 'POLYLINE', spatial_reference=spatialReference))
    arcpy.management.SplitLineAtPoint(Edge, singlePoint, saida)
    edges = []
    with arcpy.da.SearchCursor(saida, ['SHAPE@', 'Shape_Length']) as cursor:
        for row in cursor:
            edges.append(row)
    vert = Graph.Vertex(pointx, max([i for i in graph.vertices.keys()]) + 1)
    lastIndex = max([i for i in graph.edges.keys()]) + 1
    undoNewEdge = [lastIndex]
    undoNewVert = [vert.index]

    graph.vertices[vert.index] = vert
    graph.vertices[vert.index].fillEdges(str(closestID) + ';' + str(lastIndex))
    
    graph.edges[lastIndex] = Graph.Edge(edges[0][0], edges[0][1], lastIndex, graph.edges[closestID].start, vert.index)
    if str(closestID) in graph.vertices[graph.edges[closestID].start].edges:
        graph.vertices[graph.edges[closestID].start].edges.remove(str(closestID))
    if str(lastIndex) not in graph.vertices[graph.edges[closestID].start].edges:
        graph.vertices[graph.edges[closestID].start].edges.append(str(lastIndex))
    graph.edges[closestID] = Graph.Edge(edges[1][0], edges[1][1], closestID, vert.index, graph.edges[closestID].end)
    
    arcpy.management.Delete(Edge)
    arcpy.management.Delete(singlePoint)
    return graph.vertices[vert.index], undoNewEdge, undoNewVert