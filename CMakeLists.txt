cmake_minimum_required(VERSION 3.19)
project(ex2)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(ex2
        alignRand.cc
        AminoAcid.h
        Atom.cc
        Atom.h
        GeomHash.h
        HashLink.h
        HashResult.h
        macros.h
        Match.cc
        Match.h
        Matrix3.cc
        Matrix3.h
        Molecule.h
        numerics.cc
        numerics.h
        PDB.cc
        PDB.h
        RigidTrans3.cc
        RigidTrans3.h
        Rotation3.cc
        Rotation3.h
        Triangle.cc
        Triangle.h
        Vector3.cc
        Vector3.h structAlign.cpp structAlign.h)
