/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
    SpalartAllmaras turbulence model with Fv3 term and the beta
    NOTE: this is just a dummy turbulence model implementation, just to 
    create the variable nuTilda. We shouldn't call any function in this class
    Instead, we should call the correct function in DAModel/DATurbulenceModel

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasFv3FieldInversion_H
#define SpalartAllmarasFv3FieldInversion_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class SpalartAllmarasFv3FieldInversion Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasFv3FieldInversion
    : public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

    // Disallow default bitwise copy construct and assignment
    SpalartAllmarasFv3FieldInversion(const SpalartAllmarasFv3FieldInversion&);
    void operator=(const SpalartAllmarasFv3FieldInversion&);

protected:
    // Protected data

    // Fields

    /// turbulence variable for SA
    volScalarField nuTilda_;

    volVectorField UData_;

    volScalarField surfaceFriction_;

    volScalarField surfaceFrictionData_;

    volScalarField wallsquare_; 

    volScalarField volume_;  

    volScalarField pData_; 

    volScalarField USingleComponentData_;

    volScalarField Cb1_;

    volScalarField sigmaNut_;

    volScalarField Cs1_;

    volScalarField Cs2_;

    volScalarField Cw1_; 

    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;

    virtual void correctNut();

public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("SpalartAllmarasFv3FieldInversion");

    // Constructors

    //- Construct from components
    SpalartAllmarasFv3FieldInversion(
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName);

    //- Destructor
    virtual ~SpalartAllmarasFv3FieldInversion()
    {
    }

    // Member Functions

    //- Re-read model coefficients if they have changed
    virtual bool read();

    //- Return the turbulence kinetic energy
    virtual tmp<volScalarField> k() const;

    //- Return the turbulence kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const;

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "SpalartAllmarasFv3FieldInversion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
