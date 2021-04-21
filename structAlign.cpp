//
// Created by user on 21/04/2021.
//

#include "structAlign.h"
#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include "Triangle.h"
#include <iostream>
#include <fstream>
using namespace std;

void compute_trans(Match &match, GeomHash<Vector3, int> &gHash, float epsilon, Molecule<Atom> &molModel,
                   Molecule<Atom> &molTarget,
                   RigidTrans3 &rig_trans) {
    // apply rotation on each atom in the model molecule and
    // add the pairs of atoms (one from target and one from model)
    // that are close enough to the match list
    for (unsigned int i = 0; i < molModel.size(); i++) {
        Vector3 mol_atom = rig_trans * molModel[i].position(); // rotate
        // find close target molecule atoms using the hash
        HashResult<int> result;
        gHash.query(mol_atom, epsilon, result); // key is mol atom coordinate
        // check if the atoms in the result are inside the distance threshold
        // the hash is a cube shape, there can be atoms further that the threshold
        for (auto x = result.begin(); x != result.end(); x++) {
            float dist = mol_atom.dist(molTarget[*x].position());
            if (dist <= epsilon) {
                float score = (1 / (1 + dist));
                match.add(*x, i, score, score);
            }
        }
        result.clear();
    }
}


bool is_rna(Molecule<Atom> mol) {
    return mol[0].isRNABackbone();
}

int main(int argc, char *argv[]) {
    // measure the run time

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " epsilon model_pdb target_pdb" << std::endl;
        exit(1);
    }

    //********Parameters********************
    float epsilon = atof(argv[1]); // the largest allowed distance between two corresponding atoms in the alignment

    std::cout << "Epsilon: " << epsilon << std::endl;

    // read the two files into Molecule
    Molecule<Atom> molModel, molTarget;
    Molecule<Atom> molModelAll, molTargetAll;

    std::ifstream fileModelAll(argv[3]);
    std::ifstream fileTargetAll(argv[2]);

    if (!fileTargetAll) {
        std::cout << "File " << argv[3] << " does not exist." << std::endl;
        return 0;
    }
    if (!fileModelAll) {
        std::cout << "File " << argv[2] << " does not exist." << std::endl;
        return 0;
    }

    molModelAll.readPDBfile(fileModelAll);
    molTargetAll.readPDBfile(fileTargetAll);
    std::ifstream fileModel(argv[3]);
    std::ifstream fileTarget(argv[2]);
    if (!is_rna(molModelAll)) {
        molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
        molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
    } else {
        molModel.readPDBfile(fileModel, PDB::PSelector());
        molTarget.readPDBfile(fileTarget, PDB::PSelector());
    }

    // calculate center of mass
    Vector3 vectModelMass(0, 0, 0);
    for (unsigned int i = 0; i < molModel.size(); i++) {
        vectModelMass += molModel[i].position();
    }
    vectModelMass /= molModel.size();

    Vector3 vectTargetMass(0, 0, 0);
    for (unsigned int i = 0; i < molTarget.size(); i++) {
        vectTargetMass += molTarget[i].position();
    }
    vectTargetMass /= molTarget.size();

    // transform the molecules to the center of the coordinate system
    molModel += (-vectModelMass);
    molTarget += (-vectTargetMass);


    // next we insert the target molecule into hash
    // this will help us to find atoms that are close faster
    GeomHash<Vector3, int> gHash(3, epsilon); // 3 is a dimension and epsilon is the size of the hash cube
    for (unsigned int i = 0; i < molTarget.size(); i++) {
        gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
    }

    // now we try random rotations and choose the best alignment from random rotations
    unsigned int iMaxSize = 0;
    RigidTrans3 rtransBest;
    float rmsd = 0.0;
    for (int i = 0; i < molTarget.size() - 2; i++) {
        Triangle target_tr = Triangle(molTarget[i].position(), molTarget[i + 1].position(),
                                      molTarget[i + 2].position());
        cout << i << " / " <<  molTarget.size() << endl;
        for (int j = 0; j < molModel.size() - 2; j++) {

            Triangle model_tr = Triangle(molModel[j].position(), molModel[j + 1].position(),
                                         molModel[j + 2].position());
            RigidTrans3 rig_trans = target_tr | model_tr;
            // match is a class that stores the correspondence list, eg.
            // pairs of atoms, one from each molecule, that are matching
            Match match;

            compute_trans(match, gHash, epsilon, molModel, molTarget, rig_trans);
            //calculates transformation that is a little better than "rotation"
            match.calculateBestFit(molTarget, molModel);
            if (iMaxSize < match.size()) {

                iMaxSize = match.size();
                rtransBest = match.rigidTrans();
                rmsd = match.rmsd();
            }
        }
    }
    for (int i = 0; i < molModelAll.size(); i++) {
        molModelAll[i].update(rtransBest * molModelAll[i].position());
    }

    ofstream transformed("transformed.pdb");
    transformed << molModelAll;
    transformed.close();
    std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
    std::cout << "Best RMSD: " << rmsd << std::endl;
    std::cout << "Rigid Trans: " <<
              RigidTrans3(Vector3(0, 0, 0), vectTargetMass) *
              rtransBest *
              RigidTrans3(Vector3(0, 0, 0), (-vectModelMass)) << std::endl;
}

