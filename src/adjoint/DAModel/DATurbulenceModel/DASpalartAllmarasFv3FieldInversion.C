/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    This file is modified from OpenFOAM's source code
    src/TurbulenceModels/turbulenceModels/RAS/SpalartAllmaras/SpalartAllmaras.C

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

#include "DASpalartAllmarasFv3FieldInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASpalartAllmarasFv3FieldInversion, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DASpalartAllmarasFv3FieldInversion, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASpalartAllmarasFv3FieldInversion::DASpalartAllmarasFv3FieldInversion(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(modelType, mesh, daOption),
      // SA parameters
      kappa_(dimensioned<scalar>::lookupOrAddToDict(
          "kappa",
          this->coeffDict_,
          0.41)),
      Cb2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cb2",
          this->coeffDict_,
          0.622)),
      Cw2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw2",
          this->coeffDict_,
          0.3)),
      Cw3_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw3",
          this->coeffDict_,
          2.0)),
      Cv1_(dimensioned<scalar>::lookupOrAddToDict(
          "Cv1",
          this->coeffDict_,
          7.1)),
      Cv2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cv2",
          this->coeffDict_,
          5.0)),
      Cw0_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw0",
          this->coeffDict_,
          0.4)),
      // Augmented variables
      nuTilda_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("nuTilda"))),
      nuTildaRes_(
          IOobject(
              "nuTildaRes",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
#ifdef CompressibleFlow
          dimensionedScalar("nuTildaRes", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("nuTildaRes", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0),
#endif
          zeroGradientFvPatchField<scalar>::typeName),
      nuTildaResPartDeriv_(
          IOobject(
              "nuTildaResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          nuTildaRes_),
      Cb1_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("Cb1"))),
      sigmaNut_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("sigmaNut"))),      
      //Cw1_(Cb1_ / sqr(kappa_) + (1.0 + Cb2_) / sigmaNut_),
      Cw1_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("Cw1"))),
      Cs1_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("Cs1"))),
      Cs2_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("Cs2"))),
      UData_(const_cast<volVectorField&>(
          mesh.thisDb().lookupObject<volVectorField>("UData"))),
      surfaceFriction_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("surfaceFriction"))),
      surfaceFrictionData_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("surfaceFrictionData"))),
      wallsquare_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("wallsquare"))),
      volume_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("volume"))),
      pData_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("pData"))),
      USingleComponentData_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("USingleComponentData"))),
      r_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("r"))),
      fw_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("fw"))),
      y_(mesh.thisDb().lookupObject<volScalarField>("yWall"))
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// SA member functions. these functions are copied from
tmp<volScalarField> DASpalartAllmarasFv3FieldInversion::chi() const
{
    return nuTilda_ / this->nu();
}

tmp<volScalarField> DASpalartAllmarasFv3FieldInversion::fv1(
    const volScalarField& chi) const
{
    return pow( scalar(1.0) - exp(scalar(-1.0)*chi/kappa_.value()/scalar(17.0)),2);
}

tmp<volScalarField> DASpalartAllmarasFv3FieldInversion::fv2(
    const volScalarField& chi,
    const volScalarField& fv1) const
{
    return scalar(1) - chi/(scalar(1) + chi*fv1);
}

tmp<volScalarField> DASpalartAllmarasFv3FieldInversion::fv3(
    const volScalarField& chi,
    const volScalarField& fv1) const
{

    const volScalarField chiByCv2((1 / Cv2_) * chi);

    return (scalar(1) + chi * fv1)
        * (1 / Cv2_)
        * (3 * (scalar(1) + chiByCv2) + sqr(chiByCv2))
        / pow3(scalar(1) + chiByCv2);
}

tmp<volScalarField> DASpalartAllmarasFv3FieldInversion::fw(
    const volScalarField& Stilda) const
{
    //Info << "begin fw"<< nl;

    volScalarField r(
        min(
            nuTilda_
                / (max(
                       Stilda,
                       dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
                   * sqr(kappa_ * y_)),
            scalar(10.0)));
    r.boundaryFieldRef() == 0.0;

    const volScalarField Cw1(Cb1_/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_);

    const volScalarField g(r + Cw2_*(pow6(r) - r));

    const volScalarField r2(pow(r,2));

    const volScalarField F0(
        0.540*r - 0.130*r2
    );

    const volScalarField F1(
        -0.213 + 0.623*r
    );

    const volScalarField F2(
        0.049 - 0.365*r + 0.316*r2
    );
    
    scalar dC2;

    scalar cutX;
    scalar F0X;
    scalar F1X;
    scalar F2X;
    scalar dF0;
    scalar dF1;
    scalar dF2;
    scalar dFw;
    //Info << "begin fw1.0"<< nl;
    
	const volScalarField C1(Cb1_/kappa_/kappa_/Cw1);
    const volScalarField C2(1.0/sigmaNut_/Cw1);
    const volScalarField C3((1.0+ Cb2_)/sigmaNut_/Cw1);

    //Info << "begin fw1.0.1"<< nl;

    const volScalarField FwUp(C1/max(r,dimensionedScalar("SMALL", r.dimensions(), SMALL))+C2*F2/F0/F0+C3*pow(F1/F0,2));

    //Info << "begin fw1.1"<< nl;
   
    cutX = Cw0_.value();

    F0X = 0.540*cutX - 0.130*pow(cutX,2);

    F1X = -0.213 + 0.623*cutX;

    F2X = 0.049 - 0.365*cutX + 0.316*pow(cutX,2);

    //Info << "begin fw1.2"<< nl;
    
	const volScalarField FwX(C1/cutX + C2*F2X/F0X/F0X + C3*pow(F1X/F0X,2));

    //Info << "begin fw1.3"<< nl;

	// //C0 continue
	const volScalarField dC1(FwX/cutX);
	dC2 = 0.0;

    // Info << "dC1 = "<<dC1<<", dC2 = "<<dC2<<nl;
    // Info << "cutx = "<<cutX<<", FwX = "<<FwX<<nl;

    const volScalarField FwDown(dC1*r + dC2*r*r);


    volScalarField fwv(g);

    forAll(fwv,k)
    {
        
        if (r[k] <= 1.0 && r[k]>=cutX){
            fwv[k] = FwUp[k];
        }
        else if (r[k]>=0.0 && r[k]<cutX)
        {
            fwv[k] = FwDown[k];
        }
        else if (r[k]>1.0)
        {
            fwv[k] = (pow(10,2.0*Cs2_[k]-1.0)-1.0)*tanh((r[k]-1.0)*5.0/pow(10,4.0*Cs1_[k]-1.0)) + 1.0;
        }
        else if (r[k]<0.0)
        {
        }
       this->r_[k] = r[k];
	   this->fw_[k] = fwv[k];
    }

    // tmp<volScalarField> fwvv
    // (
    //     new volScalarField
    //     (   
    //         IOobject
    //         (   
    //             "fwvv",
    //             this->runTime_.constant(),
    //             this->mesh_,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),  
    //         fwv
    //     )
    // );

    return fwv * 1.0;
    //return g * pow((1.0 + pow6(Cw3_)) / (pow6(g) + pow6(Cw3_)), 1.0 / 6.0);
}

tmp<volScalarField> DASpalartAllmarasFv3FieldInversion::DnuTildaEff() const
{
    return tmp<volScalarField>(
        new volScalarField("DnuTildaEff", (nuTilda_ + this->nu()) / sigmaNut_));
}

// Augmented functions
void DASpalartAllmarasFv3FieldInversion::correctModelStates(wordList& modelStates) const
{
    /*
    Description:
        Update the name in modelStates based on the selected physical model at runtime

    Example:
        In DAStateInfo, if the modelStates reads:
        
        modelStates = {"nut"}
        
        then for the SA model, calling correctModelStates(modelStates) will give:
    
        modelStates={"nuTilda"}
        
        while calling correctModelStates(modelStates) for the SST model will give 
        
        modelStates={"k","omega"}
        
        We don't udpate the names for the radiation model becasue users are 
        supposed to set modelStates={"G"}
    */

    // replace nut with nuTilda
    forAll(modelStates, idxI)
    {
        word stateName = modelStates[idxI];
        if (stateName == "nut")
        {
            modelStates[idxI] = "nuTilda";
        }
    }
}

void DASpalartAllmarasFv3FieldInversion::correctNut()
{
    /*
    Description:
        Update nut based on other turbulence variables and update the BCs
        Also update alphat if is present
    */

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    nut_ = nuTilda_ * fv1;

    nut_.correctBoundaryConditions();

    // this is basically BasicTurbulenceModel::correctNut();
    this->correctAlphat();

    return;
}

void DASpalartAllmarasFv3FieldInversion::correctBoundaryConditions()
{
    /*
    Description:
        Update turbulence variable boundary values
    */

    // correct the BCs for the perturbed fields
    nuTilda_.correctBoundaryConditions();
}

void DASpalartAllmarasFv3FieldInversion::updateIntermediateVariables()
{
    /*
    Description:
        Update nut based on nuTilda. Note: we need to update nut and its BC since we 
        may have perturbed other turbulence vars that affect the nut values
    */

    this->correctNut();
}

void DASpalartAllmarasFv3FieldInversion::correctStateResidualModelCon(List<List<word>>& stateCon) const
{
    /*
    Description:
        Update the original variable connectivity for the adjoint state 
        residuals in stateCon. Basically, we modify/add state variables based on the
        original model variables defined in stateCon.

    Input:
    
        stateResCon: the connectivity levels for a state residual, defined in Foam::DAJacCon

    Example:
        If stateCon reads:
        stateCon=
        {
            {"U", "p", "nut"},
            {"p"}
        }
    
        For the SA turbulence model, calling this function for will get a new stateCon
        stateCon=
        {
            {"U", "p", "nuTilda"},
            {"p"}
        }
    
        For the SST turbulence model, calling this function will give
        stateCon=
        {
            {"U", "p", "k", "omega"},
            {"p", "U"}
        }
        ***NOTE***: we add a extra level of U connectivity because nut is 
        related to grad(U), k, and omega in SST!
    */

    forAll(stateCon, idxI)
    {
        forAll(stateCon[idxI], idxJ)
        {
            word conStateName = stateCon[idxI][idxJ];
            if (conStateName == "nut")
            {
                stateCon[idxI][idxJ] = "nuTilda";
            }
        }
    }
}

void DASpalartAllmarasFv3FieldInversion::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    /*
    Description:
        Add the connectivity levels for all physical model residuals to allCon

    Input:
        allCon: the connectivity levels for all state residual, defined in DAJacCon

    Example:
        If stateCon reads:
        allCon=
        {
            "URes":
            {
               {"U", "p", "nut"},
               {"p"}
            }
        }
    
        For the SA turbulence model, calling this function for will get a new stateCon,
        something like this:
        allCon=
        {
            "URes":
            {
               {"U", "p", "nuTilda"},
               {"p"}
            },
            "nuTildaRes": 
            {
                {"U", "phi", "nuTilda"},
                {"U"}
            }
        }

    */

    word pName;

    if (mesh_.thisDb().foundObject<volScalarField>("p"))
    {
        pName = "p";
    }
    else if (mesh_.thisDb().foundObject<volScalarField>("p_rgh"))
    {
        pName = "p_rgh";
    }
    else
    {
        FatalErrorIn(
            "Neither p nor p_rgh was found in mesh.thisDb()!"
            "addModelResidualCon failed to setup turbulence residuals!")
            << exit(FatalError);
    }

    // NOTE: for compressible flow, it depends on rho so we need to add T and p
#ifdef IncompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "nuTilda", "phi"}, // lv0
            {"U", "nuTilda"}, // lv1
            {"nuTilda"} // lv2
        });
#endif

#ifdef CompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "T", pName, "nuTilda", "phi"}, // lv0
            {"U", "T", pName, "nuTilda"}, // lv1
            {"T", pName, "nuTilda"} // lv2
        });
#endif
}

void DASpalartAllmarasFv3FieldInversion::correct(label printToScreen)
{
    /*
    Descroption:
        Solve the residual equations and update the state. This function will be called 
        by the DASolver. It is needed because we want to control the output frequency
        of the residual convergence every 100 steps. If using the correct from turbulence
        it will output residual every step which will be too much of information.
    */

    // We set the flag solveTurbState_ to 1 such that in the calcResiduals function
    // we will solve and update nuTilda
    solveTurbState_ = 1;
    dictionary dummyOptions;
    dummyOptions.set("printToScreen", printToScreen);
    this->calcResiduals(dummyOptions);
    // after it, we reset solveTurbState_ = 0 such that calcResiduals will not
    // update nuTilda when calling from the adjoint class, i.e., solveAdjoint from DASolver.
    solveTurbState_ = 0;
}

void DASpalartAllmarasFv3FieldInversion::calcResiduals(const dictionary& options)
{
    /*
    Descroption:
        If solveTurbState_ == 1, this function solve and update nuTilda, and 
        is the same as calling turbulence.correct(). If solveTurbState_ == 0,
        this function compute residuals for turbulence variables, e.g., nuTildaRes_

    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use the first order scheme for div(phi,nuTilda)

        p_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        nuTildaRes_: If solveTurbState_ == 0, update the residual field variable

        nuTilda_: If solveTurbState_ == 1, update nuTilda
    */

    // Copy and modify based on the "correct" function

    word divNuTildaScheme = "div(phi,nuTilda)";

    label isPC = 0;

    if (!solveTurbState_)
    {
        isPC = options.getLabel("isPC");

        if (isPC)
        {
            divNuTildaScheme = "div(pc)";
        }
    }

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    volScalarField Cw1(Cb1_/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_);

    const volScalarField Stilda(
        ::sqrt(2.0) * mag(skew(fvc::grad(U_)))
        + this->fv2(chi, fv1) * nuTilda_ / sqr(kappa_ * y_));

    tmp<fvScalarMatrix> nuTildaEqn(
        fvm::ddt(phase_, rho_, nuTilda_)
            + fvm::div(phaseRhoPhi_, nuTilda_, divNuTildaScheme)
            - fvm::laplacian(phase_ * rho_ * DnuTildaEff(), nuTilda_)
            - Cb2_ / sigmaNut_ * phase_ * rho_ * magSqr(fvc::grad(nuTilda_))
        == Cb1_ * phase_ * rho_ * Stilda * nuTilda_ 
            - fvm::Sp(Cw1 * phase_ * rho_ * fw(Stilda) * nuTilda_  / sqr(y_), nuTilda_));

    nuTildaEqn.ref().relax();

    if (solveTurbState_)
    {
        label printToScreen = options.getLabel("printToScreen");

        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverNuTilda = solve(nuTildaEqn);
        DAUtility::primalResidualControl(solverNuTilda, printToScreen, "nuTilda");

        DAUtility::boundVar(allOptions_, nuTilda_, printToScreen);
        nuTilda_.correctBoundaryConditions();

        // ***************** NOTE*****************
        // In the original SA, it is correctNut(fv1) and fv1 is not
        // updated based on the latest nuTilda. We use correctNut which
        // recompute fv1 with the latest nuTilda
        this->correctNut();
    }
    else
    {
        // calculate residuals
        nuTildaRes_ = nuTildaEqn.ref() & nuTilda_;
        // need to normalize residuals
        normalizeResiduals(nuTildaRes);
    }

    return;
}

void DASpalartAllmarasFv3FieldInversion::getTurbProdTerm(volScalarField& prodTerm) const
{
    /*
    Description:
        Return the value of the production term from the Spalart Allmaras model 
    */

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volScalarField Stilda(
         ::sqrt(2.0) * mag(skew(fvc::grad(U_)))
        + this->fv2(chi, fv1) * nuTilda_ / sqr(kappa_ * y_));
    
    volScalarField SAProdTerm = Cb1_ * phase_ * rho_ * Stilda * nuTilda_;

    forAll(SAProdTerm, cellI)
    {
        prodTerm[cellI] = SAProdTerm[cellI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
