from .mat_utils import mkvc, ndgrid, uniqueRows
import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator, interp1d
from scipy.spatial import cKDTree
import scipy.sparse as sp


def surface2ind_topo(mesh, topo, gridLoc="CC", method="nearest", fill_value=np.nan):
    """
    Get active indices from topography

    Parameters
    ----------

    :param TensorMesh mesh: TensorMesh object on which to discretize the topography
    :param numpy.ndarray topo: [X,Y,Z] topographic data
    :param str gridLoc: 'CC' or 'N'. Default is 'CC'.
                        Discretize the topography
                        on cells-center 'CC' or nodes 'N'
    :param str method: 'nearest' or 'linear' or 'cubic'. Default is 'nearest'.
                       Interpolation method for the topographic data
    :param float fill_value: default is np.nan. Filling value for extrapolation

    Returns
    -------

    :param numpy.ndarray actind: index vector for the active cells on the mesh
                               below the topography
    """
    if mesh._meshType == "TENSOR":

        if mesh.dim == 3:
            # Check if Topo points are inside of the mesh
            xmin, xmax = mesh.vectorNx.min(), mesh.vectorNx.max()
            xminTopo, xmaxTopo = topo[:, 0].min(), topo[:, 0].max()
            ymin, ymax = mesh.vectorNy.min(), mesh.vectorNy.max()
            yminTopo, ymaxTopo = topo[:, 1].min(), topo[:, 1].max()
            if (
                (xminTopo > xmin)
                or (xmaxTopo < xmax)
                or (yminTopo > ymin)
                or (ymaxTopo < ymax)
            ):
                # If not, use nearest neihbor to extrapolate them
                Ftopo = NearestNDInterpolator(topo[:, :2], topo[:, 2])
                xinds = np.logical_or(
                    xminTopo < mesh.vectorNx, xmaxTopo > mesh.vectorNx
                )
                yinds = np.logical_or(
                    yminTopo < mesh.vectorNy, ymaxTopo > mesh.vectorNy
                )
                XYOut = ndgrid(mesh.vectorNx[xinds], mesh.vectorNy[yinds])
                topoOut = Ftopo(XYOut)
                topo = np.vstack((topo, np.c_[XYOut, topoOut]))

            if gridLoc == "CC":
                XY = ndgrid(mesh.vectorCCx, mesh.vectorCCy)
                Zcc = mesh.gridCC[:, 2].reshape(
                    (np.prod(mesh.vnC[:2]), mesh.nCz), order="F"
                )
                gridTopo = griddata(
                    topo[:, :2], topo[:, 2], XY, method=method, fill_value=fill_value
                )
                actind = [gridTopo >= Zcc[:, ixy] for ixy in range(mesh.vnC[2])]
                actind = np.hstack(actind)

            elif gridLoc == "N":

                XY = ndgrid(mesh.vectorNx, mesh.vectorNy)
                gridTopo = griddata(
                    topo[:, :2], topo[:, 2], XY, method=method, fill_value=fill_value
                )
                gridTopo = gridTopo.reshape(mesh.vnN[:2], order="F")

                if mesh._meshType not in ["TENSOR", "CYL", "BASETENSOR"]:
                    raise NotImplementedError(
                        "Nodal surface2ind_topo not implemented for {0!s} mesh".format(
                            mesh._meshType
                        )
                    )

                # TODO: this will only work for tensor meshes
                Nz = mesh.vectorNz[1:]
                actind = np.array([False] * mesh.nC).reshape(mesh.vnC, order="F")

                for ii in range(mesh.nCx):
                    for jj in range(mesh.nCy):
                        actind[ii, jj, :] = [
                            np.all(gridTopo[ii : ii + 2, jj : jj + 2] >= Nz[kk])
                            for kk in range(len(Nz))
                        ]

        elif mesh.dim == 2:
            # Check if Topo points are inside of the mesh
            xmin, xmax = mesh.vectorNx.min(), mesh.vectorNx.max()
            xminTopo, xmaxTopo = topo[:, 0].min(), topo[:, 0].max()
            if (xminTopo > xmin) or (xmaxTopo < xmax):
                fill_value = "extrapolate"

            Ftopo = interp1d(topo[:, 0], topo[:, 1], fill_value=fill_value, kind=method)

            if gridLoc == "CC":
                gridTopo = Ftopo(mesh.gridCC[:, 0])
                actind = mesh.gridCC[:, 1] <= gridTopo

            elif gridLoc == "N":

                gridTopo = Ftopo(mesh.vectorNx)
                if mesh._meshType not in ["TENSOR", "CYL", "BASETENSOR"]:
                    raise NotImplementedError(
                        "Nodal surface2ind_topo not implemented for {0!s} mesh".format(
                            mesh._meshType
                        )
                    )

                # TODO: this will only work for tensor meshes
                Ny = mesh.vectorNy[1:]
                actind = np.array([False] * mesh.nC).reshape(mesh.vnC, order="F")

                for ii in range(mesh.nCx):
                    actind[ii, :] = [
                        np.all(gridTopo[ii : ii + 2] > Ny[kk]) for kk in range(len(Ny))
                    ]

        else:
            raise NotImplementedError("surface2ind_topo not implemented for 1D mesh")

    elif mesh._meshType == "TREE":
        if mesh.dim == 3:
            if gridLoc == "CC":
                # Compute unique XY location
                uniqXY = uniqueRows(mesh.gridCC[:, :2])

                if method == "nearest":
                    Ftopo = NearestNDInterpolator(topo[:, :2], topo[:, 2])
                elif method == "linear":
                    # Check if Topo points are inside of the mesh
                    xmin, xmax = mesh.x0[0], mesh.hx.sum() + mesh.x0[0]
                    xminTopo, xmaxTopo = topo[:, 0].min(), topo[:, 0].max()
                    ymin, ymax = mesh.x0[1], mesh.hy.sum() + mesh.x0[1]
                    yminTopo, ymaxTopo = topo[:, 1].min(), topo[:, 1].max()
                    if (
                        (xminTopo > xmin)
                        or (xmaxTopo < xmax)
                        or (yminTopo > ymin)
                        or (ymaxTopo < ymax)
                    ):
                        # If not, use nearest neihbor to extrapolate them
                        Ftopo = NearestNDInterpolator(topo[:, :2], topo[:, 2])
                        xinds = np.logical_or(
                            xminTopo < uniqXY[0][:, 0], xmaxTopo > uniqXY[0][:, 0]
                        )
                        yinds = np.logical_or(
                            yminTopo < uniqXY[0][:, 1], ymaxTopo > uniqXY[0][:, 1]
                        )
                        inds = np.logical_or(xinds, yinds)
                        XYOut = uniqXY[0][inds, :]
                        topoOut = Ftopo(XYOut)
                        topo = np.vstack((topo, np.c_[XYOut, topoOut]))
                    Ftopo = LinearNDInterpolator(topo[:, :2], topo[:, 2])
                else:
                    raise NotImplementedError(
                        "Only nearest and linear method are available for TREE mesh"
                    )
                actind = np.zeros(mesh.nC, dtype="bool")
                npts = uniqXY[0].shape[0]
                for i in range(npts):
                    z = Ftopo(uniqXY[0][i, :])
                    inds = uniqXY[2] == i
                    actind[inds] = mesh.gridCC[inds, 2] < z[0]
            # Need to implement
            elif gridLoc == "N":
                raise NotImplementedError("gridLoc=N is not implemented for TREE mesh")
            else:
                raise Exception("gridLoc must be either CC or N")

        elif mesh.dim == 2:

            if gridLoc == "CC":
                # Compute unique X location
                uniqX = np.unique(
                    mesh.gridCC[:, 0], return_index=True, return_inverse=True
                )

                if method == "nearest":
                    Ftopo = interp1d(topo[:, 0], topo[:, -1], kind="nearest")
                elif method == "linear":
                    # Check if Topo points are inside of the mesh
                    xmin, xmax = mesh.x0[0], mesh.hx.sum() + mesh.x0[0]
                    xminTopo, xmaxTopo = topo[:, 0].min(), topo[:, 0].max()
                    if (xminTopo > xmin) or (xmaxTopo < xmax):
                        # If not, use nearest neihbor to extrapolate them
                        Ftopo = interp1d(topo[:, 0], topo[:, -1], kind="nearest")
                        xinds = np.logical_or(
                            xminTopo < uniqX[0][:, 0], xmaxTopo > uniqX[0][:, 0]
                        )
                        XOut = uniqX[0][xinds, :]
                        topoOut = Ftopo(XOut)
                        topo = np.vstack((topo, np.c_[XOut, topoOut]))
                    Ftopo = interp1d(topo[:, 0], topo[:, -1], kind="nearest")
                else:
                    raise NotImplementedError(
                        "Only nearest and linear method are available for TREE mesh"
                    )
                actind = np.zeros(mesh.nC, dtype="bool")
                npts = uniqX[0].shape[0]
                for i in range(npts):
                    z = Ftopo(uniqX[0][i])
                    inds = uniqX[2] == i
                    actind[inds] = mesh.gridCC[inds, 1] < z
            # Need to implement
            elif gridLoc == "N":
                raise NotImplementedError("gridLoc=N is not implemented for TREE mesh")
            else:
                raise Exception("gridLoc must be either CC or N")

        else:
            raise NotImplementedError("surface2ind_topo not implemented for 1D mesh")

    return mkvc(actind)


def surface_layer_index(mesh, topo, index=0):
    """
        Find the ith layer below topo
    """

    actv = np.zeros(mesh.nC, dtype="bool")
    # Get cdkTree to find top layer
    tree = cKDTree(mesh.gridCC)

    def ismember(a, b):
        bind = {}
        for i, elt in enumerate(b):
            if elt not in bind:
                bind[elt] = i
        return np.vstack([bind.get(itm, None) for itm in a])

    grid_x, grid_y = np.meshgrid(mesh.vectorCCx, mesh.vectorCCy)
    zInterp = mkvc(
        griddata(topo[:, :2], topo[:, 2], (grid_x, grid_y), method="nearest")
    )

    # Get nearest cells
    r, inds = tree.query(np.c_[mkvc(grid_x), mkvc(grid_y), zInterp])
    inds = np.unique(inds)

    # Extract vertical neighbors from Gradz operator
    Dz = mesh._cellGradzStencil
    Iz, Jz, _ = sp.find(Dz)
    jz = np.sort(Jz[np.argsort(Iz)].reshape((int(Iz.shape[0] / 2), 2)), axis=1)
    for ii in range(index):

        members = ismember(inds, jz[:, 1])
        inds = np.squeeze(jz[members, 0])

    actv[inds] = True

    return actv
