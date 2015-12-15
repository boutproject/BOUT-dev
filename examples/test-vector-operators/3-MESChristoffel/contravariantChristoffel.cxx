// ******** Checking the Christoffel symbols **************
#include <bout/physicsmodel.hxx>

class ChristoffelCheck : public PhysicsModel {
private:
    // The Christoffel symbols
    Field3D G1_11;
    Field3D G1_12;
    Field3D G1_13;
    Field3D G1_22;
    Field3D G1_23;
    Field3D G1_33;
    Field3D G2_11;
    Field3D G2_12;
    Field3D G2_13;
    Field3D G2_22;
    Field3D G2_23;
    Field3D G2_33;
    Field3D G3_11;
    Field3D G3_12;
    Field3D G3_13;
    Field3D G3_22;
    Field3D G3_23;
    Field3D G3_33;

protected:
    int init(bool restarting) {
        // Save the Christoffel symbols
        G1_11 = mesh->G1_11;
        G1_12 = mesh->G1_12;
        G1_13 = mesh->G1_13;
        G1_22 = mesh->G1_22;
        G1_23 = mesh->G1_23;
        G1_33 = mesh->G1_33;
        G2_11 = mesh->G2_11;
        G2_12 = mesh->G2_12;
        G2_13 = mesh->G2_13;
        G2_22 = mesh->G2_22;
        G2_23 = mesh->G2_23;
        G2_33 = mesh->G2_33;
        G3_11 = mesh->G3_11;
        G3_12 = mesh->G3_12;
        G3_13 = mesh->G3_13;
        G3_22 = mesh->G3_22;
        G3_23 = mesh->G3_23;
        G3_33 = mesh->G3_33;
        SAVE_ONCE6(G1_11, G1_12, G1_13, G1_22, G1_23, G1_33);
        SAVE_ONCE6(G2_11, G2_12, G2_13, G2_22, G2_23, G2_33);
        SAVE_ONCE6(G3_11, G3_12, G3_13, G3_22, G3_23, G3_33);
        dump.write();
        dump.close();
        output << "\nWriting complete!" << std::endl;
        // Wait for all processors to write data
        MPI_Barrier(BoutComm::get());
        return 0;
    }

    int rhs(BoutReal t) {
        // Doesn't do anything
        return 0;
    }
};

BOUTMAIN(ChristoffelCheck);
