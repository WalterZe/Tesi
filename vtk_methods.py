#!/usr/bin/env python

import sys,math,os
import vtk
from vmtk import vtkvmtk
from vmtk import pypes

#METODI GESTIONE FILE
def ReadPolyData(filename):
    """ Legge un file contenete PolyData"""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()
    
def ReadSTL(filename):
    """ Legge un file contente una gemetria con estensione .STL"""
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()
def ReadTxt(filename):
    """ Legge dal file i punti e li restituisce come lista
    """
    points =[]
    p = [0, 0, 0]
    with open(filename) as file:
        linea = file.readline()
        for linea in file:
            p_str = linea.split()
            for i in range(3):
                p[i] = float(p_str[i])
            points.append([p[0], p[1], p[2]])
    return points
    
def WritePolyData(surface, filename):
    """ Scrive PolyData su file """
    if not os.path.exists(filename):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInput(surface)
        writer.SetFileName(filename)
        writer.Write()
        return filename
    else:
        return filename
    
def WritePointsTxt(points, fn):
    with open(fn, "w") as file:
        for i in range(len(points)): 
            string=''
            for k in range(3): 
                string+=str(points[i][k])+' '
            string+='\n'
            print string
            file.write(string)
        file.close()
    return None
        
#METODI GESTIONE RENDERER
def CreateActor(oggetto):
    """ Crea un attore """
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(oggetto)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor

def CreateSphere(center, radius=0.1, color=(0.0, 0.0, 0.0)): 
    """crea una sfera color r,g,b"""
    if center==None:
        print 'ERRORE: vtk_methods.CreateSphere: parametri mancanti'
        pass
    else:
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(center)
        sphere.SetRadius(radius)
        sphere.SetThetaResolution(18)
        sphere.SetPhiResolution(18)
        actor_sphere = CreateActor(sphere.GetOutput())
        prop_shere = actor_sphere.GetProperty()
        prop_shere.SetColor(color)
        return actor_sphere

def CreateArrow(center, color=(0.0, 0.0, 0.0)):
    if center==None:
        print 'ERRORE: vtk_methods.CreateArrow: parametri mancanti'
        pass
    else:        
        vertGlyph = vtk.vtkGlyph3D()
        vertGlyph.SetSource(center)
        vertGlyph.ScalingOn()
        vertGlyph.SetScaleModeToScaleByVector()
        vertGlyph.SetScaleFactor(0.6)
        vertGlyph.OrientOn()
        vertGlyph.SetVectorModeToUseVector()
        glyph = vtkGlyph3D()
        glyph.SetInputConnection(vertex_geom.GetOutputPort())
        glyph.SetSourceConnection(0, cube.GetOutputPort())
        actor_glyph = CreateActor(vertGlyph.GetOutput())
        prop_glyph = actor_glyph.GetProperty()
        prop_glyph.SetColor(color)
        return actor_glyph
    
def CreateRenderMulti(actor, actor1, actor2, actor3):
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(actor1)
    renderer.AddActor(actor2)
    renderer.AddActor(actor3)
    return renderer
    
def CreateRenderWindow(renderer, x, y):
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetPosition(x, y)
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.LightFollowCameraOff()
    interactorStyle = vtk.vtkInteractoStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactorStyle)
    renderWindow.SetInteractor(renderWindowInteractor)
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()
    
def RenderWindowMulti(ren):
    renderWindow = vtk.vtkRenderWindow()
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.LightFollowCameraOff()
    interactorStyle = vtk.vtkInteractoStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactorStyle)
    renderWindow.SetInteractor(renderWindowInteractor)
   # Definisco i range dei viewport
    xmins=[0, .5, 0, .5]
    xmaxs=[0.5, 1, 0.5, 1]
    ymins=[0, 0, .5, .5]
    ymaxs=[0.5, 0.5, 1, 1]
    for i in range(4):
        renderWindow.AddRenderer(ren)
        ren.SetViewport(xmins[i], ymins[i], xmaxs[i], ymaxs[i])
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindow.SetWindowName('RW: Mutiple ViewPorts')
    renderWindowInteractor.Start()
    
def CreateSurface(input, midPoint):
    connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
    connectivityFilter.SetInput(input)
    connectivityFilter.SetClosestPoint(midPoint)
    connectivityFilter.SetExtractionModeToClosestPointRegion()
    connectivityFilter.Update()
    contour = connectivityFilter.GetOutput()
    numberOfContourPoints = contour.GetNumberOfPoints()
    surf = vtk.vtkPolyData()
    sectionPoints = vtk.vtkPoints()
    sectionCellArray = vtk.vtkCellArray() 
    sectionCellArray.InsertNextCell(numberOfContourPoints)
    for j in range(numberOfContourPoints):
        point = contour.GetPoint(j)
        sectionPoints.InsertNextPoint(point)
        sectionCellArray.InsertCellPoint(j)
    sectionCellArray.InsertCellPoint(0)
    surf.SetPoints(sectionPoints)
    surf.SetPolys(sectionCellArray)
    return surf

def ComputeCenterlinePype(fn_in, fn_out):
    str0='vmtkcenterlines -ifile '+fn_in+' -endpoints 1 -seedselector carotidprofiles'
    str0= 'vmtkcenterlines -ifile '+fn_in+' -endpoints 1'
    str1=' -costfunction 1/R^2 -resampling 1 -resamplingstep 2.5 -ofile '+fn_out
    myArguments=str0+str1
    myPype=pypes.PypeRun(myArguments)
    return fn_out

def CreateCoords_versore(o, r):
    """ Ritorna una lista di attori contenenti i il sistema di coordinate:
    o = origine
    r = versore"""
    points = []
    Lines=[]
    Polygon = vtk.vtkPolyData()
    ac=[]
    
    points = vtk.vtkPoints()
            
    points.SetNumberOfPoints(4)

    points.SetPoint(0, self.midPoint)
    points.SetPoint(1, [self.FrenetBinormalArray[0]+self.midPoint[0], self.FrenetBinormalArray[1]+self.midPoint[1], self.FrenetBinormalArray[2]+self.midPoint[2]])
    points.SetPoint(2, [self.FrenetNormalArray[0]+self.midPoint[0], self.FrenetNormalArray[1]+self.midPoint[1], self.FrenetNormalArray[2]+self.midPoint[2]])
    points.SetPoint(3, [self.FrenetTangentArray[0]+self.midPoint[0], self.FrenetTangentArray[1]+self.midPoint[1], self.FrenetTangentArray[2]+self.midPoint[2]])
     
    points.SetPoint(0, o)
    points.SetPoint(1, [o[0]+r[0], o[1]         , o[2]])
    points.SetPoint(2, [o[0]        , o[1]+r[1] ,  o[2]])
    points.SetPoint(3, [o[0]        , o[1]         ,  o[2]+r[2]])
    
    polyLine0 = vtk.vtkPolyLine()
    polyLine0.GetPointIds().SetNumberOfIds(2)
    polyLine0.GetPointIds().SetId(0,0)
    polyLine0.GetPointIds().SetId(1,1)
        
    polyLine1 = vtk.vtkPolyLine()
    polyLine1.GetPointIds().SetNumberOfIds(2)
    polyLine1.GetPointIds().SetId(0,0)
    polyLine1.GetPointIds().SetId(1,2)
    
    polyLine2 = vtk.vtkPolyLine()
    polyLine2.GetPointIds().SetNumberOfIds(2)
    polyLine2.GetPointIds().SetId(0,0)
    polyLine2.GetPointIds().SetId(1,3)
    
    cells0 = vtk.vtkCellArray()
    cells0.InsertNextCell(polyLine0)
    cells0.InsertNextCell(polyLine1)
    cells0.InsertNextCell(polyLine2)
    
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells0)
    
    ac=[]
        
    
    ac.append(CreateSphere(points.GetPoint(0), 0.05, [1, 1, 1]))
    ac.append(CreateSphere(points.GetPoint(1), 0.1, [1, 0, 0]))
    ac.append(CreateSphere(points.GetPoint(2), 0.1, [0, 1, 0]))
    ac.append(CreateSphere(points.GetPoint(3), 0.1, [0, 0, 1]))
        
                
    ac.append(CreateActor(polyData))
    return ac
    
def CreateCoords(o, x, y, z):
        """ Ritorna una lista di attori di un sistema di coordinate cartesiane """
        points = []
        Lines=[]
        Polygon = vtk.vtkPolyData()
        ac=[]
        
        points = vtk.vtkPoints()
                
        points.SetNumberOfPoints(4)

#        points.SetPoint(0, self.midPoint)
#        points.SetPoint(1, [self.FrenetBinormalArray[0]+self.midPoint[0], self.FrenetBinormalArray[1]+self.midPoint[1], self.FrenetBinormalArray[2]+self.midPoint[2]])
#        points.SetPoint(2, [self.FrenetNormalArray[0]+self.midPoint[0], self.FrenetNormalArray[1]+self.midPoint[1], self.FrenetNormalArray[2]+self.midPoint[2]])
#        points.SetPoint(3, [self.FrenetTangentArray[0]+self.midPoint[0], self.FrenetTangentArray[1]+self.midPoint[1], self.FrenetTangentArray[2]+self.midPoint[2]])
         
        points.SetPoint(0, o)
        points.SetPoint(1, [o[0]+x[0], o[1]+x[1],  o[2]+x[2]])
        points.SetPoint(2, [o[0]+y[0], o[1]+y[1],  o[2]+y[2]])
        points.SetPoint(3, [o[0]+z[0], o[1]+z[1],  o[2]+z[2]])
        
        polyLine0 = vtk.vtkPolyLine()
        polyLine0.GetPointIds().SetNumberOfIds(2)
        polyLine0.GetPointIds().SetId(0,0)
        polyLine0.GetPointIds().SetId(1,1)
            
        polyLine1 = vtk.vtkPolyLine()
        polyLine1.GetPointIds().SetNumberOfIds(2)
        polyLine1.GetPointIds().SetId(0,0)
        polyLine1.GetPointIds().SetId(1,2)
        
        polyLine2 = vtk.vtkPolyLine()
        polyLine2.GetPointIds().SetNumberOfIds(2)
        polyLine2.GetPointIds().SetId(0,0)
        polyLine2.GetPointIds().SetId(1,3)
        
        cells0 = vtk.vtkCellArray()
        cells0.InsertNextCell(polyLine0)
        cells0.InsertNextCell(polyLine1)
        cells0.InsertNextCell(polyLine2)
        
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cells0)
        
        ac=[]
            
        
        ac.append(CreateSphere(points.GetPoint(0), 0.05, [1, 1, 1]))
        ac.append(CreateSphere(points.GetPoint(1), 0.1, [1, 0, 0]))
        ac.append(CreateSphere(points.GetPoint(2), 0.1, [0, 1, 0]))
        ac.append(CreateSphere(points.GetPoint(3), 0.1, [0, 0, 1]))
            
                    
        ac.append(CreateActor(polyData))
        return ac


def CreateVersor(o, r, color=[0, 0, 0]):
    """ Ritorna una lista di attori contenetni un versore
   o = origine
  r= versore """
    points = []
    Lines=[]
    Polygon = vtk.vtkPolyData()
    ac=[]
    
    points = vtk.vtkPoints()
            
    points.SetNumberOfPoints(2)
    
    r = normalize(r)
         
    points.SetPoint(0, o)
    points.SetPoint(1, [o[0]+r[0], o[1]+r[1], o[2]+r[2]])
    
    polyLine0 = vtk.vtkPolyLine()
    polyLine0.GetPointIds().SetNumberOfIds(2)
    polyLine0.GetPointIds().SetId(0,0)
    polyLine0.GetPointIds().SetId(1,1)
    
    cells0 = vtk.vtkCellArray()
    cells0.InsertNextCell(polyLine0)
    
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells0)
    
    ac=[]
        
    
    ac.append(CreateSphere(points.GetPoint(0), 0.05, [1, 1, 1]))
    ac.append(CreateSphere(points.GetPoint(1), 0.1, color))
    
    ac.append(CreateActor(polyData))
    return ac
    
#METODI GEOMETRICI

#vettori
def add(x, y): return [a+b for a, b in zip(x, y)]
def diff(x, y): return [a-b for a, b in zip(x, y)]
def mul(t, v): return [t*a for a in v]

def prodottovettoriale(v, w):
    (v1, v2, v3) = v; (w1, w2, w3) = w
    return [v2*w3-v3*w2,  v3*w1-v1*w3,  v1*w2-v2*w1]
    
def prodottoscalare (u, v):
    s = 0
    for x, y in zip(u, v): s+=x*y
    return s


def modulo(a):
    """ calcola il modulo di un vettore """
    return math.sqrt(prodottoscalare(a, a))

def norm(u, v):
    """ calcola la norma di due vettori """
    s=0
    for a, b in zip(u, v) : s+=(a-b)**2 
    return math.sqrt(s)

def min(a, b):
    if a<b:
        return a
    else:
        return b



                     

def centerOfMass(a):
    """ calcola il centro di masssa di una lista di punti"""
    meanPoint = [0, 0, 0]
    for i in len(a):
        meanPoint[0]= meanPoint[0] + a[0]
        meanPoint[1]= meanPoint[1] + a[1]
        meanPoint[2]= meanPoint[2] + a[2]
    meanPoint[0]=meanPoint[0]/len(a)
    meanPoint[1]=meanPoint[1]/len(a)
    meanPoint[2]=meanPoint[2]/len(a)
    return meanPoint

def somma_vettori(a, b):
    """ somma i vettori """
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
    
def dot(b, a):
    return [a[0]*b, a[1]*b, a[2]*b]
def vettore(a, b):
    """ crea un vettore tra due punti """
    return [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
    
def versors(a):
    """ ritorna i versori di un vettore """
    return [a[0]/modulo(a), a[1]/modulo(a), a[2]/modulo(a)]
def versore(a, b):
    return normalize(vettore(a, b))
def rettaSpazio2Punti(P1, P2, t):
    """ fornisce l'equazione dei punti x,y,z della retta passante per il punto P e parallela al vettore v
    """
    v = [P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2] ]
    
    return (v[0]*t+P1[0], v[1]*t+P1[1], v[2]*t+P1[2])
def proiezionePuntosuPiano(P, n, p): 
    c = prodottoscalare(P, n)
    v = mul((c-prodottoscalare(n, p))/float(prodottoscalare(n, n)), n) 
    return add(p, v)

    
    
# QUATERNIONI
def normalize ( v,  tolerance = 0.00001):
    mag2 = prodottoscalare(v, v)
    if abs(mag2 -1.0) > tolerance:
        mag = math.sqrt(mag2)
        v = [v[0]/mag, v[1]/mag, v[2]/mag]
    return v
    
def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 -  x1* x2  -  y1 * y2 -  z1* z2
    x =  w1 * x2 + x1* w2 + y1 * z2 -  z1* y2
    y =  w1 * y2 + y1* w1 + z1 * x2 -  x1* z2
    z =  w1 * z2 + z1 * w2 + x1 * y2 - y1* x2
    return w, x, y, z
    
def q_conj(q):
    q = normalize(q)
    w, x, y, z = q
    return(w, -x, -y, -z)
    
def qv_mult(q1, v1):
    v1 = normalize(v1)
    q2 = (0.0, ) + v1
    return q_mult(q_mult(q1, q2),  q_conj(q1))[1:]
    
def axisangleToQ(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = cos(theta)
    x = x* sin(theta)
    y = y*sin(theta)
    z = z*sin(theta)
    return w, x, y, z
    
def QtoAxisAngle(q):
    w, v = q[0],  q[1:]
    theta = acos(w)*2.0
    return normalize(v), theta
    
# Generatori di geometrie
def creaCerchio(radius, n):
    """ Crea un cerchio sul piano x,y con raggio e numero di punti dati
    ritorna la lista dei punti
    """
    # equazione cerchio:
    theta = (2*math.pi)/n
    x = lambda r, theta: round(r*math.cos(theta), 4)
    y = lambda r, theta: round(r*math.sin(theta), 4)
    points=[]
    for i in range(n+1):
        points.append([x(radius, theta*i) ,  y(radius, theta*i) , 0])
    return points
    
WritePointsTxt(creaCerchio(2, 50), '/home/walter/Desktop/cerchiotest.txt')
