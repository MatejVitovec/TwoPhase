#include <cmath>
#include <iostream>

#include "Limiter.hpp"

/*Field<Vars<5>> Limiter::calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    //Field<Vars<5>> out(cells.size());
    Field<Vars<5>> out(cells.size(), Vars<5>({10.0, 10.0, 10.0, 10.0, 10.0}));

    Field<Vars<5>> wCmax(cells.size());
    Field<Vars<5>> wCmin(cells.size());

    //over faces
    for (int i = 0; i < faces.size(); i++)
    {
        int ownerIndex = owners[i];  

        wCmax[ownerIndex] = max(wCmax[ownerIndex], wr[i]);
        wCmin[ownerIndex] = min(wCmin[ownerIndex], wr[i]);

        int neighbourIndex = neighbours[i];
        if (neighbourIndex >= 0)
        {
            wCmax[neighbourIndex] = max(wCmax[neighbourIndex], wl[i]);
            wCmin[neighbourIndex] = min(wCmin[neighbourIndex], wl[i]);
        }
    }

    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbours[i];

        Compressible wCOwner = wl[i];
        Vars<5> denominatorOwner = dot(grad[owner], faces[i].midpoint - cells[owner].center);
        Vars<5> phiCnOwner;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorOwner[k] > 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmax[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else if (denominatorOwner[k] < 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmin[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else
            {
                phiCnOwner[k] = 1.0;
            }                
        }

        out[owner] = min(out[owner], phiCnOwner);

        if (neighbor < 0)
        {
            continue;
        }
        
        Compressible wCNeighbor = wr[i];
        Vars<5> denominatorNeighbor = dot(grad[neighbor], faces[i].midpoint - cells[neighbor].center);
        Vars<5> phiCnNeighbor;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorNeighbor[k] > 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmax[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else if (denominatorNeighbor[k] < 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmin[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else
            {
                phiCnNeighbor[k] = 1.0;
            }                
        }

        out[neighbor] = min(out[neighbor], phiCnNeighbor);
    }
    return out;
}*/

/*Field<Vars<5>> Limiter::calculateLimiter(const Field<Compressible>& w, const std::vector<std::vector<Compressible>>& boundaryField, Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    //Field<Vars<5>> out(cells.size());
    Field<Vars<5>> out(cells.size(), Vars<5>({10.0, 10.0, 10.0, 10.0, 10.0}));

    Field<Vars<5>> wCmax(cells.size());
    Field<Vars<5>> wCmin(cells.size());

    //over faces
    for (int i = 0; i < faces.size(); i++)
    {
        int ownerIndex = owners[i];
        int neighbourIndex = neighbours[i];

        //wCmax[ownerIndex] = max(wCmax[ownerIndex], wr[i]);
        //wCmin[ownerIndex] = min(wCmin[ownerIndex], wr[i]);

        //if (neighbourIndex >= 0)
        //{
        //    wCmax[neighbourIndex] = max(wCmax[neighbourIndex], wl[i]);
        //    wCmin[neighbourIndex] = min(wCmin[neighbourIndex], wl[i]);
        //}

        //////////
        if (neighbourIndex >= 0)
        {
            wCmax[ownerIndex] = max(wCmax[ownerIndex], w[neighbourIndex]);
            wCmin[ownerIndex] = min(wCmin[ownerIndex], w[neighbourIndex]);
            wCmax[neighbourIndex] = max(wCmax[neighbourIndex], w[ownerIndex]);
            wCmin[neighbourIndex] = min(wCmin[neighbourIndex], w[ownerIndex]);
        }
    }

    for (int boundaryId = 0; boundaryId < boundaryField.size(); boundaryId++)
    {
        const std::vector<int>& boundaryFaces = mesh.getBoundaryList()[boundaryId].facesIndex;
        for (int i = 0; i < boundaryFaces.size(); i++)
        {
            wCmax[owners[boundaryFaces[i]]] = max(wCmax[owners[boundaryFaces[i]]], boundaryField[boundaryId][i]);
            wCmin[owners[boundaryFaces[i]]] = min(wCmin[owners[boundaryFaces[i]]], boundaryField[boundaryId][i]);
        }
    }
    
    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbours[i];

        //Compressible wCOwner = wl[i];
        Compressible wCOwner = w[owner];
        Vars<5> denominatorOwner = dot(grad[owner], faces[i].midpoint - cells[owner].center);
        Vars<5> phiCnOwner;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorOwner[k] > 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmax[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else if (denominatorOwner[k] < 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmin[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else
            {
                phiCnOwner[k] = 1.0;
            }                
        }

        out[owner] = min(out[owner], phiCnOwner);

        if (neighbor < 0)
        {
            continue;
        }
        
        //Compressible wCNeighbor = wr[i];
        Compressible wCNeighbor = w[neighbor];
        Vars<5> denominatorNeighbor = dot(grad[neighbor], faces[i].midpoint - cells[neighbor].center);
        Vars<5> phiCnNeighbor;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorNeighbor[k] > 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmax[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else if (denominatorNeighbor[k] < 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmin[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else
            {
                phiCnNeighbor[k] = 1.0;
            }                
        }

        out[neighbor] = min(out[neighbor], phiCnNeighbor);
    }

    return out;
}*/

Field<Vars<5>> Limiter::calculateLimiter(const VolField<Compressible>& w, Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    //Field<Vars<5>> out(cells.size());
    Field<Vars<5>> out(cells.size(), Vars<5>({10.0, 10.0, 10.0, 10.0, 10.0}));

    Field<Vars<5>> wCmax(cells.size());
    Field<Vars<5>> wCmin(cells.size());

    //over faces
    for (int i = 0; i < faces.size(); i++)
    {
        int ownerIndex = owners[i];
        int neighbourIndex = neighbours[i];

        /*wCmax[ownerIndex] = max(wCmax[ownerIndex], wr[i]);
        wCmin[ownerIndex] = min(wCmin[ownerIndex], wr[i]);

        if (neighbourIndex >= 0)
        {
            wCmax[neighbourIndex] = max(wCmax[neighbourIndex], wl[i]);
            wCmin[neighbourIndex] = min(wCmin[neighbourIndex], wl[i]);
        }*/

        //////////
        if (neighbourIndex >= 0)
        {
            wCmax[ownerIndex] = max(wCmax[ownerIndex], w[neighbourIndex]);
            wCmin[ownerIndex] = min(wCmin[ownerIndex], w[neighbourIndex]);
            wCmax[neighbourIndex] = max(wCmax[neighbourIndex], w[ownerIndex]);
            wCmin[neighbourIndex] = min(wCmin[neighbourIndex], w[ownerIndex]);
        }
    }

    for (int boundaryId = 0; boundaryId < w.boundarySize(); boundaryId++)
    {
        const std::vector<int>& boundaryFaces = mesh.getBoundaryList()[boundaryId].facesIndex;
        for (int i = 0; i < boundaryFaces.size(); i++)
        {
            wCmax[owners[boundaryFaces[i]]] = max(wCmax[owners[boundaryFaces[i]]], w.boundary(boundaryId)[i]);
            wCmin[owners[boundaryFaces[i]]] = min(wCmin[owners[boundaryFaces[i]]], w.boundary(boundaryId)[i]);
        }
    }
    
    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbours[i];

        //Compressible wCOwner = wl[i];
        Compressible wCOwner = w[owner];
        Vars<5> denominatorOwner = dot(grad[owner], faces[i].midpoint - cells[owner].center);
        Vars<5> phiCnOwner;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorOwner[k] > 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmax[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else if (denominatorOwner[k] < 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmin[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else
            {
                phiCnOwner[k] = 1.0;
            }                
        }

        out[owner] = min(out[owner], phiCnOwner);

        if (neighbor < 0)
        {
            continue;
        }
        
        //Compressible wCNeighbor = wr[i];
        Compressible wCNeighbor = w[neighbor];
        Vars<5> denominatorNeighbor = dot(grad[neighbor], faces[i].midpoint - cells[neighbor].center);
        Vars<5> phiCnNeighbor;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorNeighbor[k] > 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmax[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else if (denominatorNeighbor[k] < 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmin[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else
            {
                phiCnNeighbor[k] = 1.0;
            }                
        }

        out[neighbor] = min(out[neighbor], phiCnNeighbor);
    }

    return out;
}



Field<Vars<9>> Limiter::calculateLimiter(const VolField<CompressibleMixture>& w, Field<Mat<9,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Vars<9>> out(cells.size(), Vars<9>({10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}));

    Field<Vars<9>> wCmax(cells.size());
    Field<Vars<9>> wCmin(cells.size());

    //over faces
    for (int i = 0; i < faces.size(); i++)
    {
        int ownerIndex = owners[i];
        int neighbourIndex = neighbours[i];

        if (neighbourIndex >= 0)
        {
            wCmax[ownerIndex] = max(wCmax[ownerIndex], w[neighbourIndex]);
            wCmin[ownerIndex] = min(wCmin[ownerIndex], w[neighbourIndex]);
            wCmax[neighbourIndex] = max(wCmax[neighbourIndex], w[ownerIndex]);
            wCmin[neighbourIndex] = min(wCmin[neighbourIndex], w[ownerIndex]);
        }
    }

    for (int boundaryId = 0; boundaryId < w.boundarySize(); boundaryId++)
    {
        const std::vector<int>& boundaryFaces = mesh.getBoundaryList()[boundaryId].facesIndex;
        for (int i = 0; i < boundaryFaces.size(); i++)
        {
            wCmax[owners[boundaryFaces[i]]] = max(wCmax[owners[boundaryFaces[i]]], w.boundary(boundaryId)[i]);
            wCmin[owners[boundaryFaces[i]]] = min(wCmin[owners[boundaryFaces[i]]], w.boundary(boundaryId)[i]);
        }
    }
    
    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbours[i];

        CompressibleMixture wCOwner = w[owner];
        Vars<9> denominatorOwner = dot(grad[owner], faces[i].midpoint - cells[owner].center);
        Vars<9> phiCnOwner;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorOwner[k] > 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmax[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else if (denominatorOwner[k] < 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmin[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else
            {
                phiCnOwner[k] = 1.0;
            }                
        }

        out[owner] = min(out[owner], phiCnOwner);

        if (neighbor < 0)
        {
            continue;
        }
        
        //Compressible wCNeighbor = wr[i];
        CompressibleMixture wCNeighbor = w[neighbor];
        Vars<9> denominatorNeighbor = dot(grad[neighbor], faces[i].midpoint - cells[neighbor].center);
        Vars<9> phiCnNeighbor;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorNeighbor[k] > 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmax[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else if (denominatorNeighbor[k] < 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmin[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else
            {
                phiCnNeighbor[k] = 1.0;
            }                
        }

        out[neighbor] = min(out[neighbor], phiCnNeighbor);
    }

    return out;
}        



double Limiter::limiterFunction(double y) const
{
    return 0.0;
}