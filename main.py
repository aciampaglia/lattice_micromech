from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np
import os
import random

# Helper function to check defect intereference
def intersects_existing_spheres(defects, x0, y0, z0, d0):
    print("{}".format(d0))
    x, y, z, d = np.array(defects['x0']), np.array(defects['y0']), np.array(defects['z0']), np.array(defects['d'])
    for i in range(len(defects['x0'])):
        distance_squared = (x[i] - x0) ** 2 + (y[i] - y0) ** 2 + (z[i] - z0) ** 2
        if distance_squared < (d0/2 + d[i]/2) ** 2:
            return True
    return False

# Helper function to create a spherical part
def create_sphere(d):
    # Create the defect
    temp_sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    temp_sketch.ConstructionLine(point1= (0.0, -10*d), point2=(0.0, 10*d))
    # Create the semi-circle for circumference
    temp_sketch.ArcByCenterEnds(center=( 0.0, 0.0), direction=CLOCKWISE, point1=(0.0, d/2), point2=(0.0, -d/2))
    temp_sketch.Line(point1=(0.0, d/2), point2=(0.0, -d/2))
    # Revolve the generator by 360
    defect_part = model.Part(dimensionality=THREE_D, name='defect', type=DEFORMABLE_BODY)
    defect_part.BaseSolidRevolve(angle=360.0, 
        flipRevolveDirection=OFF, sketch= temp_sketch)
    del temp_sketch
    print('Sphere created')
    return defect_part

# Helper function to cut the defect from the RVE and update the defects dict
def generate_defect(defect_part, position):
    global rveInstance, defects
    x0, y0, z0 = position
    print("Creating defect in in: ({:.1f}, {:.1f}, {:.1f})".format(x0, y0, z0))
    defectInstance = model.rootAssembly.Instance(dependent=ON, name=
        'defectInstance', part=defect_part)
    # Translate the defect in x0, y0, z0
    model.rootAssembly.translate(instanceList=('defectInstance',),vector=(x0, y0, z0))
    # Cut the defect from the RVE part and update it
    rveInstance =   model.rootAssembly.InstanceFromBooleanCut(
                    cuttingInstances = (defectInstance, ),
                    instanceToBeCut = rveInstance,
                    name='RVE', 
                    originalInstances=DELETE)
    defects['d'].append(d); defects['x0'].append(x0); defects['y0'].append(y0); defects['z0'].append(z0);
    print("Defect created")

# Generate a spherical (and periodic) defect in the RVE at random (or specified) location
def add_defect(d, L, position='random', periodic=True):
    global rveInstance, defects
    # Define the position
    if position == 'random':
        x0, y0, z0 = random.uniform(-L/2, L/2), random.uniform(-L/2, L/2), random.uniform(-L/2, L/2)
        while intersects_existing_spheres(defects, x0, y0, z0, d):
            x0, y0, z0 = random.uniform(-L/2, L/2), random.uniform(-L/2, L/2), random.uniform(-L/2, L/2)
        position = [x0, y0, z0]
    # Generate the spherical periodic void and locate in the position
    defect_part = create_sphere(d)
    generate_defect(defect_part, position)
    # new_part = model.parts['RVE']
    # Cut the periodic defect is periodic=True
    if periodic:    #Generate a periodic defect in the RVE
        if (L/2-abs(x0))<d/2 and (L/2-abs(y0))>d/2 and (L/2-abs(z0))>d/2:               #sfere su facce normali asse X
            print("Defect interesecting the x face")
            periodic_position = [-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))), y0, z0]
            defect_part = create_sphere(d)
            generate_defect(defect_part, periodic_position)
        elif (L/2-abs(y0))<d/2 and (L/2-abs(x0))>d/2 and (L/2-abs(z0))>d/2:             #sfere su facce normali asse Y
            print("Defect interesecting the y face")
            periodic_position = [x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))), z0]
            defect_part = create_sphere(d)
            generate_defect(defect_part, periodic_position)
        elif (L/2-abs(z0))<d/2 and (L/2-abs(x0))>d/2 and (L/2-abs(y0))>d/2:             #sfere su facce normali asse Z
            print("Defect interesecting the z face")
            periodic_position = [x0, y0, -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))]
            defect_part = create_sphere(d)
            generate_defect(defect_part, periodic_position)
        elif (L/2-abs(z0))>d/2 and (L/2-abs(x0))<d/2 and (L/2-abs(y0))<d/2:             #sfere su bordo tra facce normali a x e y
            print("Defect interesecting the x-y edge")
            ry=L/2-abs(y0)
            rx=L/2-abs(x0)
            rsq=sqrt(rx**2+ry**2)
            if rsq > d/2:
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))),y0 , z0])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , z0])
            else:
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))), -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , z0])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))),y0 , z0])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , z0])
        elif (L/2-abs(z0))<d/2 and (L/2-abs(x0))<d/2 and (L/2-abs(y0))>d/2:    #sfere su bordo tra facce normali a x e z
            print("Defect interesecting the x-z edge")
            rz=L/2-abs(z0)
            rx=L/2-abs(x0)
            rsq=sqrt(rx**2+rz**2)
            if rsq > d/2:
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))),y0 , z0])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
            else:
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))), y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))),y0 , z0])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
        elif (L/2-abs(z0))<d/2 and (L/2-abs(x0))>d/2 and (L/2-abs(y0))<d/2:    #sfere su bordo tra facce normali a y e z
            print("Defect interesecting the y-z edge")
            ry=L/2-abs(y0)
            rz=L/2-abs(z0)
            rsq=sqrt(rz**2+ry**2)
            if rsq > d/2:
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0 ,-np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))),z0])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
            else:
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0,y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
                defect_part = create_sphere(d)
                generate_defect(defect_part, position=[x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , z0])
        elif (L/2-abs(z0))<d/2 and (L/2-abs(x0))<d/2 and (L/2-abs(y0))<d/2:    #sfere sui vertici del RVE
            print("Defect interesecting the corners")
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[x0,y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[x0, -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , z0])
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))), y0 , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))),y0 , z0])
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))), -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , z0])
            defect_part = create_sphere(d)
            generate_defect(defect_part, position=[-np.sign(x0)*(abs(x0)+2*(L/2-abs(x0))), -np.sign(y0)*(abs(y0)+2*(L/2-abs(y0))) , -np.sign(z0)*(abs(z0)+2*(L/2-abs(z0)))])
    return position

# Helper function to refien the mesh around a defect
def create_mesh(global_size = 0.05, refined_size=0.01):
    global defects
    for i in range(len(defects['x0'])):
        x0, y0, z0 = defects['x0'][i], defects['y0'][i], defects['z0'][i]
        d = defects['d'][i]
        # Create partitions for better meshing
        xyPlane = model.parts['RVE'].DatumPlaneByPrincipalPlane(offset= z0, principalPlane=XYPLANE)
        yzPlane = model.parts['RVE'].DatumPlaneByPrincipalPlane(offset= x0, principalPlane=YZPLANE)
        xzPlane = model.parts['RVE'].DatumPlaneByPrincipalPlane(offset= y0, principalPlane=XZPLANE)
        for i in range(2, len(model.parts['RVE'].datums)+2):
            model.parts['RVE'].PartitionCellByDatumPlane(cells= model.parts['RVE'].cells, datumPlane= model.parts['RVE'].datums[i])
        model.parts['RVE'].setMeshControls(elemShape=TET, regions=
            model.parts['RVE'].cells, technique=FREE)
        model.parts['RVE'].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.05)
    
        r45 = d/2*np.sin(np.deg2rad(45))
        for radial_distance in [[r45, r45], [r45, -r45], [-r45, r45], [-r45, -r45]]:
            model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
                deviationFactor=0.1, edges=
                model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0+radial_distance[1], z0),)), minSizeFactor=0.1, size=0.01)
            model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
                deviationFactor=0.1, edges=
                model.parts['RVE'].edges.findAt(((x0, y0+radial_distance[0], z0+radial_distance[1]),)), minSizeFactor=0.1, size=0.01)
            model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
                deviationFactor=0.1, edges=
                model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0, z0+radial_distance[1]),)), minSizeFactor=0.1, size=0.01)

# Generate a cubic RVE with spherical defects
def generate_cubic_rve(d, L, position='random', target_vf=False):
    # Create the cubic rve
    temp_sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    temp_sketch.rectangle(point1=(-L/2,-L/2),point2=(L/2,L/2))
    rve_part = model.Part(dimensionality=THREE_D, name='RVE', type=DEFORMABLE_BODY)
    rve_part.BaseSolidExtrude(sketch=temp_sketch, depth=L)
    del temp_sketch
    rveInstance = model.rootAssembly.Instance(dependent=ON, name='rve-1', part=rve_part)
    model.rootAssembly.translate(instanceList=('rve-1',), vector=(0.0, 0.0, -L/2))
    vf = 1.0
    # Add the first defect
    defects = {'d': [], 'x0':[], 'y0':[], 'z0':[]}
    defect_location = add_defect(d, L, position)
    vf -= 4/3*np.pi*(d/2)**3/L**3
    # Add more defect if there is a target volume fraction
    if target_vf:
        while vf>target_vf:
            defect_location = add_defect(d, L, position='random')
            vf -= 4/3*np.pi*(d/2)**3/L**3
    # Create the mesh
    r45 = d/2*np.sin(np.deg2rad(45))
    for radial_distance in [[r45, r45], [r45, -r45], [-r45, r45], [-r45, -r45]]:
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0+radial_distance[1], z0),)), minSizeFactor=0.1, size=0.01)
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0, y0+radial_distance[0], z0+radial_distance[1]),)), minSizeFactor=0.1, size=0.01)
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0, z0+radial_distance[1]),)), minSizeFactor=0.1, size=0.01)
    rve_part.generateMesh()
    rve_part.setElementType(elemTypes=(ElemType(
        elemCode=C3D8R, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
        elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)), regions=(rve_part.cells.getSequenceFromMask((
        '[#ff ]', ), ), ))
    # Create an empty set to store matching nodes
    tolerance = 1e-3
    allNodes = new_part.nodes
    node_region = allNodes.getByBoundingBox(-L, L/2-tolerance, -L, L, L/2+tolerance, L) 
    new_part.Set(name="Top", nodes=node_region)
    tolerance = 1e-3
    allNodes = new_part.nodes
    node_region = allNodes.getByBoundingBox(-L, -L/2-tolerance, -L, L, -L/2+tolerance, L) 
    new_part.Set(name="Bottom", nodes=node_region)
    model.rootAssembly.regenerate()

def generate_cylindrical_rve(d, D, L, x0, y0, z0):
    # Create the sketch
    temp_sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    temp_sketch.ConstructionLine(point1= (0.0, -L), point2=(0.0, L))
    # Create the semi-circle for circumference
    temp_sketch.ArcByCenterEnds(center=( 0.0, 0.0), direction=CLOCKWISE, point1=(0.0, d/2), point2=(0.0, -d/2))
    temp_sketch.Line(point1=(0.0, d/2), point2=(0.0, -d/2))
    # Revolve the generator by 360
    defect_part = model.Part(dimensionality=THREE_D, name='defect', type=DEFORMABLE_BODY)
    defect_part.BaseSolidRevolve(angle=360.0, 
        flipRevolveDirection=OFF, sketch= temp_sketch)
    del temp_sketch
    # Create the RVE
    temp_sketch = model.ConstrainedSketch(name='__profile__', sheetSize= 200.0)
    temp_sketch.ConstructionLine(point1= (0.0, -L), point2=(0.0, L))
    temp_sketch.rectangle(point1=(0.0, -L/2), point2=(D/2, L/2))
    rve_part = model.Part(dimensionality=THREE_D, name='rve', type= DEFORMABLE_BODY)
    rve_part.BaseSolidRevolve(angle=360.0, 
        flipRevolveDirection=OFF, sketch= temp_sketch)
    del temp_sketch
    # Create the assembly
    defectInstance = model.rootAssembly.Instance(dependent=ON, name=
        'defect-1', part=defect_part)
    rveInstance = model.rootAssembly.Instance(dependent=ON, name='rve-1'
        , part=rve_part)
    # Localize the defect in x0, y0, z0
    model.rootAssembly.translate(instanceList=('defect-1', 
        ), vector=(x0, y0, z0))
    model.rootAssembly.InstanceFromBooleanCut(
        cuttingInstances = (defectInstance, ),
        instanceToBeCut = rveInstance,
        name='RVE', 
        originalInstances=DELETE)
    new_part = model.parts['RVE']
    # Create partitions for better meshing
    xyPlane = new_part.DatumPlaneByPrincipalPlane(offset= z0, principalPlane=XYPLANE)
    yzPlane = new_part.DatumPlaneByPrincipalPlane(offset= x0, principalPlane=YZPLANE)
    xzPlane = new_part.DatumPlaneByPrincipalPlane(offset= y0, principalPlane=XZPLANE)
    print("There are {} datum planes".format(len(new_part.datums)))
    for i in range(2, len(new_part.datums)+2):
        new_part.PartitionCellByDatumPlane(cells= new_part.cells, datumPlane= new_part.datums[i])
    new_part.setMeshControls(elemShape=TET, regions=
        new_part.cells, technique=FREE)
    new_part.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.01)
    r45 = d/2*np.sin(np.deg2rad(45))
    for radial_distance in [[r45, r45], [r45, -r45], [-r45, r45], [-r45, -r45]]:
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0+radial_distance[1], z0),)), minSizeFactor=0.1, size=0.005)
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0, y0+radial_distance[0], z0+radial_distance[1]),)), minSizeFactor=0.1, size=0.005)
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0, z0+radial_distance[1]),)), minSizeFactor=0.1, size=0.005)
    for distance in [d/2+0.1, -d/2-0.1]:
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0+distance, y0, z0),)), minSizeFactor=0.1, size=refined_size)
#       model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
#           deviationFactor=0.1, edges=
#           model.parts['RVE'].edges.findAt(((x0, y0+distance, z0),)), minSizeFactor=0.1, size=refined_size)
        model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=
            model.parts['RVE'].edges.findAt(((x0, y0, z0+distance),)), minSizeFactor=0.1, size=refined_size)
    new_part.generateMesh()
    new_part.setElementType(elemTypes=(ElemType(
        elemCode=C3D8R, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
        elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)), regions=(new_part.cells.getSequenceFromMask((
        '[#ff ]', ), ), ))
    # Create an empty set to store matching nodes
    tolerance = 1e-3
    allNodes = new_part.nodes
    node_region = allNodes.getByBoundingBox(-D, L/2-tolerance, -D, D, L/2+tolerance, D) 
    new_part.Set(name="Top", nodes=node_region)
    tolerance = 1e-3
    allNodes = new_part.nodes
    node_region = allNodes.getByBoundingBox(-D, -L/2-tolerance, -D, D, -L/2+tolerance, D) 
    new_part.Set(name="Bottom", nodes=node_region)
    model.rootAssembly.regenerate()

def create_dogbone_shape(d, Dg, Dtot, Lg, Ltot, R, x0, y0, z0, global_size=0.1, refined_size=0.05):
	# Create the sketch
	temp_sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
	temp_sketch.ConstructionLine(point1= (0.0, -Ltot), point2=(0.0, Ltot))
	# Create the semi-circle for circumference
	temp_sketch.ArcByCenterEnds(center=( 0.0, 0.0), direction=CLOCKWISE, point1=(0.0, d/2), point2=(0.0, -d/2))
	temp_sketch.Line(point1=(0.0, d/2), point2=(0.0, -d/2))
	# Revolve the generator by 360
	defect_part = model.Part(dimensionality=THREE_D, name='defect', type=DEFORMABLE_BODY)
	defect_part.BaseSolidRevolve(angle=360.0, 
		flipRevolveDirection=OFF, sketch= temp_sketch)
	del temp_sketch
	# Create the RVE
	s = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
	s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
	s.Line(point1=(Dg/2, -Lg/2), point2=(Dg/2, Lg/2))
	y2 = Lg/2+np.sqrt(R**2-(R+Lg/2-Ltot/2)**2)
	s.ArcByStartEndTangent(point1=(Dg/2, Lg/2), point2=(Dtot/2, y2), vector=(0., 1.))
	s.ArcByStartEndTangent(point1=(Dg/2, -Lg/2), point2=(Dtot/2, -y2), vector=(0., -1.))
	s.Line(point1=(Dtot/2, y2), point2=(Dtot/2, Ltot/2))
	s.Line(point1=(Dtot/2, Ltot/2), point2=(0.0, Ltot/2))
	s.Line(point1=(0.0, Ltot/2), point2=(0.0, -Ltot/2))
	s.Line(point1=(0.0, -Ltot/2), point2=(Dtot/2, -Ltot/2))
	s.Line(point1=(Dtot/2, -Ltot/2), point2=(Dtot/2, -y2))
	# Revolve the generator by 360
	rve_part = model.Part(dimensionality=THREE_D, name='rve', type=DEFORMABLE_BODY)
	rve_part.BaseSolidRevolve(angle=360.0, 
		flipRevolveDirection=OFF, sketch= s)
	del s
	# Create the assembly
	defectInstance = model.rootAssembly.Instance(dependent=ON, name=
		'defect-1', part=defect_part)
	rveInstance = model.rootAssembly.Instance(dependent=ON, name='rve-1'
		, part=rve_part)
	# Localize the defect in x0, y0, z0
	model.rootAssembly.translate(instanceList=('defect-1', 
		), vector=(x0, y0, z0))
	model.rootAssembly.InstanceFromBooleanCut(
		cuttingInstances = (defectInstance, ),
		instanceToBeCut = rveInstance,
		name='RVE', 
		originalInstances=DELETE)
	new_part = model.parts['RVE']
	# Create partitions for better meshing
	xyPlane = new_part.DatumPlaneByPrincipalPlane(offset= z0, principalPlane=XYPLANE)
	yzPlane = new_part.DatumPlaneByPrincipalPlane(offset= x0, principalPlane=YZPLANE)
	xzPlane = new_part.DatumPlaneByPrincipalPlane(offset= y0, principalPlane=XZPLANE)
	print("There are {} datum planes".format(len(new_part.datums)))
	for i in range(2, len(new_part.datums)+2):
		new_part.PartitionCellByDatumPlane(cells= new_part.cells, datumPlane= new_part.datums[i])
	new_part.setMeshControls(elemShape=TET, regions=
		new_part.cells, technique=FREE)
	new_part.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=global_size)
	r45 = d/2*np.sin(np.deg2rad(45))
	# Refine sphere edges
	for radial_distance in [[r45, r45], [r45, -r45], [-r45, r45], [-r45, -r45]]:
		model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
			deviationFactor=0.1, edges=
			model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0+radial_distance[1], z0),)), minSizeFactor=0.1, size=refined_size)
		model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
			deviationFactor=0.1, edges=
			model.parts['RVE'].edges.findAt(((x0, y0+radial_distance[0], z0+radial_distance[1]),)), minSizeFactor=0.1, size=refined_size)
		model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
			deviationFactor=0.1, edges=
			model.parts['RVE'].edges.findAt(((x0+radial_distance[0], y0, z0+radial_distance[1]),)), minSizeFactor=0.1, size=refined_size)
	for distance in [d/2+0.1, -d/2-0.1]:
		model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
			deviationFactor=0.1, edges=
			model.parts['RVE'].edges.findAt(((x0+distance, y0, z0),)), minSizeFactor=0.1, size=refined_size)
#		model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
#			deviationFactor=0.1, edges=
#			model.parts['RVE'].edges.findAt(((x0, y0+distance, z0),)), minSizeFactor=0.1, size=refined_size)
		model.parts['RVE'].seedEdgeBySize(constraint=FINER, 
			deviationFactor=0.1, edges=
			model.parts['RVE'].edges.findAt(((x0, y0, z0+distance),)), minSizeFactor=0.1, size=refined_size)
	new_part.generateMesh()
	new_part.setElementType(elemTypes=(ElemType(
		elemCode=C3D8R, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
		elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD, 
		secondOrderAccuracy=OFF, distortionControl=DEFAULT)), regions=(new_part.cells.getSequenceFromMask((
		'[#ff ]', ), ), ))
	# Create an empty set to store matching nodes
	tolerance = 0.001
	allNodes = new_part.nodes
	node_region = allNodes.getByBoundingBox(-Dtot-tolerance, Ltot/2-tolerance, -Dtot-tolerance, Dtot+tolerance, Ltot/2+tolerance, Dtot+tolerance) 
	new_part.Set(name="Top", nodes=node_region)
	allNodes = new_part.nodes
	node_region = allNodes.getByBoundingBox(-Dtot-tolerance, -Ltot/2-tolerance, -Dtot-tolerance, Dtot+tolerance, -Ltot/2+tolerance, Dtot+tolerance) 
	new_part.Set(name="Bottom", nodes=node_region)
	model.rootAssembly.regenerate()


def assign_material(mat_name, rho, E, v):
    mat = model.Material(name=mat_name)
    mat.Density(table=((rho, ), ))
    mat.Elastic(table=((E, v), ))
    mat.Plastic(table=(
    (207.1, 0),
    (225.0,   0.000602076),
    (237.0, 0.001183677),
    (246.0, 0.001781526),
    (256.8, 0.002608351),
    (273.9, 0.004945053),
    (287.1, 0.007262566),
    (306.1, 0.010563982),
    (318.3, 0.012795855),
    (344.3, 0.018513212),
    (358.9, 0.02358372),
    (376.5, 0.030505689),
    (392.6, 0.039538212),
    (409.3, 0.056435717),
    (411.7, 0.066352881),
    (415.1, 0.097022066)))
    sec = model.HomogeneousSolidSection(name='Solid', material=mat_name, thickness=None)
    model.parts['RVE'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(cells=model.parts['RVE'].cells),
        sectionName='Solid', thicknessAssignment=FROM_SECTION)

def apply_bcs(uy):
# Generate step and output field
    model.StaticStep(name='Static_1', previous='Initial')
    model.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'E', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'EVOL'))
    model.fieldOutputRequests['F-Output-1'].setValues(numIntervals=35)
    
    # Contrain bottom part
    model.DisplacementBC(amplitude=UNSET, createStepName=
        'Static_1', distributionType=UNIFORM, fieldName='', fixed=OFF, 
        localCsys= None, name='BC-1', 
        region= model.rootAssembly.instances['RVE-1'].sets['Bottom'], 
        u1= UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=UNSET, ur3=SET)    
    
    # Apply displacement on top
    model.DisplacementBC(amplitude=UNSET, createStepName=
        'Static_1', distributionType=UNIFORM, fieldName='', fixed=OFF, 
        localCsys= None, name='BC-2', 
        region= model.rootAssembly.instances['RVE-1'].sets['Top'], 
        u1= UNSET, u2=uy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

def run_job(model_name, ncpu):
    job = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model=model_name, modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=model_name, nodalOutputPrecision=
        SINGLE, numCpus=ncpu, numDomains=ncpu, numGPUs=0, queue=None, resultsFormat=ODB, 
        scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    job.writeInput()
    job.submit(consistencyChecking=OFF)
    return job


## ## Test single RVE
## model_name='prova_2'
## model = mdb.Model(modelType=STANDARD_EXPLICIT, name=model_name)
## d, D, L = 0.35, 1.1, 2
## radial_location = np.random.rand(1)*(D/2-d)*0.85 #radial location random
## x0, y0, z0 = radial_location[0], 0.0, 0.0
## generate_cylindrical_rve(d, D, L, x0, y0, z0)
## mat_name = 'alsi10mg'
## rho, E, v = 2.7, 53e3, 0.33
## assign_material(mat_name, rho, E, v)
## apply_bcs(0.2)
## job = run_job(model_name, 8)

## CY-RVE size analysis
output_directory="C:\Users\Administrator\OneDrive - Politecnico di Torino\14_additive_lattice\simulazioni_rve\cyrve_dmax_3x3"
if not os.path.isdir(output_directory):
    os.mkdir(output_directory)
os.chdir(output_directory)
# The RVE size convergence is repeated for the 5 largest defects with a random radial location, fixed in the mid-section
rve_size = [1.4] # mm
defects_file = 'dmax_defects_2x2.txt'
defect_values = np.loadtxt(defects_file, skiprows=1, usecols=3)*1e-3
repetitions = 1
for i, d in enumerate(defect_values):
	print("Generating defect with size: {}".format(d))
	if not os.path.isdir('defect_{}'.format(int(d*100))):
		os.mkdir('defect_{}'.format(int(d*100)))
	os.chdir('defect_{}'.format(int(d*100)))
	for size in rve_size:
		for rep in range(repetitions):
			model_name='size_analysis_{}_{}_{}'.format(int(d*100), int(size*100), rep+1)
			if not os.path.isdir(model_name):
				os.mkdir(model_name)
			os.chdir(model_name)
			model = mdb.Model(modelType=STANDARD_EXPLICIT, name=model_name)
			D, L = size, 1.0
			radial_location = 0.0 #np.random.rand(1)*(D/2-d)*0.85 #radial location random
			x0, y0, z0 = 0.0, 0.0, 0.0
			generate_cylindrical_rve(d, D, L, x0, y0, z0)
			mat_name = 'alsi10mg'
			rho, E, v = 2.7, 53e3, 0.33
			assign_material(mat_name, rho, E, v)
			apply_bcs(0.2)
			job = run_job(model_name, 16)
			job.waitForCompletion()
			os.chdir('..')
	os.chdir('..')
