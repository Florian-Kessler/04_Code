import numpy as np
import pyvista as pv
import SimpleITK as sitk
import skimage.morphology as morph
from plyer import notification
import sys


def stl2Mask(dict_, path, filename, resolution, mask_name='BoneTest.mhd',
             controlplot=False, reshape=True, closing=True):
    """
    This function creates a mask form any stl file and returns a 3d array mask - and store the mask as mhd in the given
    path.
    :param dict_: empty dictionary to store all mask information
    :param path: path to store a mhd file of the mask
    :param filename: name of the stl file
    :param resolution: resolution (voxel size) of the mask
    :param mask_name: name of the mhd file
    :param controlplot: If true a control 3d image of the stl will pop up - close it to proceed
    :param reshape: sometimes the order of the slices for the 3d array does not match - activate it if mhd looks wierd
    :param closing: sometimes there are some small holes in the mask - activate it if needed
    :return: 3d array mask
    """

    # read in the stl to generate the mask
    reader = pv.get_reader(filename)
    mesh = reader.read()

    if controlplot:
        mesh.plot()
        voxels = pv.voxelize(mesh, density=resolution)
        p = pv.Plotter()
        p.add_mesh(voxels, color=True, show_edges=True)
        p.add_mesh(mesh, color="red", opacity=0.5)
        p.show()

    x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
    x = np.arange(x_min, x_max, resolution)
    y = np.arange(y_min, y_max, resolution)
    z = np.arange(z_min, z_max, resolution)
    x, y, z = np.meshgrid(x, y, z)

    # Create unstructured grid from the structured grid
    grid = pv.StructuredGrid(x, y, z)
    ugrid = pv.UnstructuredGrid(grid)

    # get part of the mesh within the mesh's bounding surface.
    selection = ugrid.select_enclosed_points(mesh.extract_surface(), tolerance=0.0, check_surface=False)
    mask_ = selection.point_data['SelectedPoints'].view(np.bool_)

    if reshape:
        # sometimes the order of the matrix gets changed
        mask = mask_.reshape([z.shape[2], z.shape[1], z.shape[0]])
        mask = mask[:, ::-1, :]
        mask = np.rot90(mask, k=-1, axes=(1, 2))
    else:
        mask = mask_.reshape(x.shape)

    mask = np.array(mask)

    if closing:
        # close small holes and gabs:
        for i in range(0, mask.shape[0]):
            mask[i, :, :] = morph.closing(mask[i, :, :], np.ones([3, 3]))
            # mask_copy[i, :, :] = morph.dilation(mask_copy[i, :, :], np.ones([3, 3]))
            # mask_copy[i, :, :] = morph.erosion(mask_copy[i, :, :], np.ones([2, 2]))

    origin = [-5.75, -8.75, -45]
    spacing = np.array([1, 1, 1]) * resolution

    mask_trans = mask.astype(np.short)
    itkmask = sitk.GetImageFromArray(mask_trans, isVector=None)
    itkmask.SetSpacing(spacing)
    itkmask.SetOrigin(origin)

    path_to_local_folder = path
    sitk.WriteImage(itkmask, f'{path_to_local_folder}/{mask_name}')

    # set bone values
    dict_['MASK_array'] = mask_trans.T
    # Needed to move the COS to the right position inside the mask
    dict_['MaskX'] = np.array([x_min, x_max])
    dict_['MaskY'] = np.array([y_min, y_max])
    dict_['MaskZ'] = np.array([z_min, z_max])

    print('BoneMeshMask')
    return dict_


def main():
    path2stl = '/home/biomech/Documents/01_Icotec/01_Experiments/99_Others/Screw_Ostp_new_20230605'
    stl = 'LowStiffnessScrew_FreeCAD_v0.stl'

    # a = vtk.vtkSTLReader()
    # a.SetFileName(f'{path2stl}/{stl}')
    # a.Update()
    # a = a.GetOutput()

    imp_info = {}
    imp_info = stl2Mask(imp_info, path2stl, f'{path2stl}/{stl}', 0.0606995, mask_name='LowStiffnessScrew_test_v0.mhd')
    notification.notify(
        title="Finished executing " + sys.argv[0],
        message="Successful",
    )
    return imp_info


if __name__ == '__main__':
    main()
