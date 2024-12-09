#ifndef CELL_HPP
#define CELL_HPP

#include <vector>
#include <memory>
#include "Face.hpp"
#include "../Vars.hpp"

class Cell
{
    public:
        enum cellType{GENERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID};

        Cell() : type(GENERAL) {};
        Cell(std::vector<int> nodesIdx) : type(GENERAL) {};
        Cell(std::vector<int> nodesIdx, cellType cType) : nodesIndex(nodesIdx), type(cType) {};

        void update(const std::vector<Face>& faceList);

        std::vector<Face> createFaces(); //for GMSH

        int getVtkType() const;

        virtual ~Cell();

        std::vector<int> nodesIndex;
        std::vector<int> ownFaceIndex;
        std::vector<int> neighborFaceIndex;

        Vars<3> center;
        Vars<3> projectedArea;
        double volume;

        void calcApproxCenter(const std::vector<Vars<3>>& nodeList);
        
    protected:
        int type;

};

inline std::ostream& operator<<(std::ostream& os, const Cell& cell)
{
    os << cell.nodesIndex.size() << " ";
    for (auto & nodeIndex : cell.nodesIndex)
    {
        os << nodeIndex << " ";
    }
    
    return os;
}


#endif // CELL_HPP