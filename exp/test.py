import os
import sys
import vtk
root =  os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


iren = vtk.vtkRenderWindowInteractor()
iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

renWin = vtk.vtkRenderWindow()
renWin.SetSize(1000, 1000)
iren.SetRenderWindow(renWin)

ren = vtk.vtkRenderer()
ren.SetBackground(0, 0, 0)
renWin.AddRenderer(ren)


def MakeActor(polydata):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)


    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor

def addPointOrder(polydata):
    nPoints = polydata.GetNumberOfPoints()

    scalars = vtk.vtkIntArray()
    scalars.SetNumberOfComponents(1)
    scalars.SetNumberOfTuples(nPoints)

    for i in range(nPoints):
        scalars.SetTuple(i, [i])

    polydata.GetPointData().SetScalars(scalars)

if __name__ == "__main__":


    
    app_path = os.path.join(root, "build", "manifold")
    input_path = "./sampleData/mandibular_nonmanifold.obj"
    output_path = "output.obj"
    resolution = 100
    
    cmd = "%s %s %s %i" % (app_path, input_path, output_path, resolution)

    
    #Write manifold obj file
    os.system(cmd)


    #input
    objReader=  vtk.vtkOBJReader()
    objReader.SetFileName(input_path)
    objReader.Update()
    inputPolyData = objReader.GetOutput()
    inputActor = MakeActor(inputPolyData)
    ren.AddActor(inputActor)


    #output
    objReader=  vtk.vtkOBJReader()
    objReader.SetFileName(output_path)
    objReader.Update()
    outputPolyData = objReader.GetOutput()
    outputActor = MakeActor(outputPolyData)
    outputActor.SetPosition(inputActor.GetBounds()[1] - inputActor.GetBounds()[0], 0, 0)
    ren.AddActor(outputActor)

    # addPointOrder(inputPolyData)
    # inputActor.GetMapper().SetScalarRange(0, inputPolyData.GetNumberOfPoints())
    
    # addPointOrder(outputPolyData)
    # outputActor.GetMapper().SetScalarRange(0, outputPolyData.GetNumberOfPoints())

    print("input : ", inputPolyData.GetNumberOfPoints(), inputPolyData.GetNumberOfCells())
    print("output :", outputPolyData.GetNumberOfPoints(), outputPolyData.GetNumberOfCells())
    print("resolution : ", resolution)

    renWin.Render()

    iren.Initialize()
    iren.Start()