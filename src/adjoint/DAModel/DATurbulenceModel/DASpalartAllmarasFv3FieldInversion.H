/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
    Child class for the SpalartAllmarasFv3 model with a betaFieldInversion field multiplying 
    to the production term. This betaFieldInversion term can then be trained to improve the
    SA turbulence model.

    This file is modified from OpenFOAM's source code
    src/TurbulenceModels/turbulenceModels/RAS/SpalartAllmaras/SpalartAllmaras.H

    OpenFOAM: The Open Source CFD Toolbox

    Copyright (C): 2011-2016 OpenFOAM Foundation

    OpenFOAM License:

        OpenFOAM is free software: you can redistribute it and/or modify it
        under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
    
        OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
        ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
        FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
        for more details.
    
        You should have received a copy of the GNU General Public License
        along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef DASpalartAllmarasFv3FieldInversion_H
#define DASpalartAllmarasFv3FieldInversion_H

#include "DATurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASpalartAllmarasFv3FieldInversion Declaration
\*---------------------------------------------------------------------------*/

class DASpalartAllmarasFv3FieldInversion
    : public DATurbulenceModel
{

protected:
    /// \name SA parameters
    //@{
    dimensionedScalar kappa_;
    dimensionedScalar Cb2_;
    dimensionedScalar Cw2_;
    dimensionedScalar Cw3_;
    dimensionedScalar Cv1_;
    dimensionedScalar Cv2_;
    dimensionedScalar Cw0_;
    //@}

    /// \name SA functions
    //@{
    tmp<volScalarField> chi() const;

    tmp<volScalarField> fv1(const volScalarField& chi) const;

    tmp<volScalarField> fv2(
        const volScalarField& chi,
        const volScalarField& fv1) const;

    tmp<volScalarField> fv3(
        const volScalarField& chi,
        const volScalarField& fv1) const;

    tmp<volScalarField> fw(const volScalarField& Stilda) const;

    //@}

    /// \name Augmented variables for residual computation
    //@{
    volScalarField& nuTilda_;
    volScalarField nuTildaRes_;
    volScalarField nuTildaResPartDeriv_;

    volScalarField r_;
    volScalarField fw_;
    //@}

    /// A coefficient Field
    volScalarField& Cb1_;

    /// A coefficient Field
    volScalarField& sigmaNut_;

    /// A coefficient Field
    volScalarField& Cw1_;

    /// A coefficient Field
    volScalarField& Cs1_; 

    /// A coefficient Field
    volScalarField& Cs2_;    

    /// reference velocity field
    volVectorField UData_; 

    /// a surface friction 'field' when using skin friction data for field inversion
    volScalarField surfaceFriction_;

    /// the reference field for surfaceFriction
    volScalarField surfaceFrictionData_;

    /// wall square
    volScalarField wallsquare_;

    /// volume
    volScalarField volume_;  

    /// the reference field for pressure
    volScalarField pData_; 

    /// reference field for profile data
    volScalarField USingleComponentData_; 

    /// 3D wall distance
    const volScalarField& y_;

    /// whether to solve for turb states
    label solveTurbState_ = 0;

    

public:
    TypeName("SpalartAllmarasFv3FieldInversion");
    // Constructors

    //- Construct from components
    DASpalartAllmarasFv3FieldInversion(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption);

    //- Destructor
    virtual ~DASpalartAllmarasFv3FieldInversion()
    {
    }

    // Member functions

    /// Return the effective diffusivity for nuTilda
    tmp<volScalarField> DnuTildaEff() const;

    //包淳新添加
    //tmp<volScalarField> f1() const;
    //tmp<volScalarField> f2() const;

    /// update the turbulence state for DAStateInfo::regStates_
    virtual void correctModelStates(wordList& modelStates) const;

    /// update nut based on other turbulence variables and update the BCs
    virtual void correctNut();

    /// update turbulence variable boundary values
    virtual void correctBoundaryConditions();

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the original variable connectivity for the adjoint state residuals in stateCon
    virtual void correctStateResidualModelCon(List<List<word>>& stateCon) const;

    /// add the model residual connectivity to stateCon
    virtual void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;

    /// compute the turbulence residuals
    virtual void calcResiduals(const dictionary& options);

    /// solve the residual equations and update the state
    virtual void correct(label printToScreen);

    /// return the value of the production term from the turbulence model 
    virtual void getTurbProdTerm(volScalarField& prodTerm) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
