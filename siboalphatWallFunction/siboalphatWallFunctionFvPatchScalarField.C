/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "siboalphatWallFunctionFvPatchScalarField.H"
#include "compressibleTurbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

siboalphatWallFunctionFvPatchScalarField::siboalphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
	Ks_(0.0)//added
{}


siboalphatWallFunctionFvPatchScalarField::siboalphatWallFunctionFvPatchScalarField
(
    const siboalphatWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
	Ks_(ptf.Ks_)//added
{}


siboalphatWallFunctionFvPatchScalarField::siboalphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85)),
	Ks_(dict.lookupOrDefault<scalar>("Ks", 0.0))//added
{}


siboalphatWallFunctionFvPatchScalarField::siboalphatWallFunctionFvPatchScalarField
(
    const siboalphatWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
	Ks_(awfpsf.Ks_)//added
{}


siboalphatWallFunctionFvPatchScalarField::siboalphatWallFunctionFvPatchScalarField
(
    const siboalphatWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
	Ks_(awfpsf.Ks_)//added
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void siboalphatWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    /*const compressibleTurbulenceModel& turbModel =
        db().lookupObject<compressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                compressibleTurbulenceModel::propertiesName,
                internalField().group()
            )
        );*/

	const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<scalarField> tnutw = turbModel.nut(patchi);
	const scalarField& nutw = tnutw();
	
	//*added
	//const scalarField& nuw = turbModel.nu().boundaryField()[patchi];
	//const tmp<scalarField> nuw = turbModel.nu(patchi);
	const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();	
	//const scalarField& alphaw = turbModel.alpha().boundaryField()[patchi];
	//const tmp<scalarField> alphaw = turbModel.alpha(patchi);
	const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    const scalarField& alphaw = talphaw();
	const fvPatchVectorField& Up = turbModel.U().boundaryField()[patchi];
	const scalarField magUp = mag(Up.patchInternalField() - Up);
	scalarField magFaceGradU = mag(Up.snGrad());
	const scalar smallNum = 1e-14;
	//tmp<scalarField> talpha(new scalarField(*this));
	//scalarField& alphaTurb = talpha();
	scalarField& alphaTurb = *this;//not sure
	scalar aFun = 0.5;
	//*

	forAll(alphaTurb, facei)
	{
		scalar Pr = muw[facei]/alphaw[facei];
		scalar tau = (muw[facei]+rhow[facei]*nutw[facei])*magFaceGradU[facei];
		scalar uTau = sqrt(max(tau/rhow[facei],smallNum));
		scalar KsPlus = uTau*Ks_/(muw[facei]/rhow[facei]);
		scalar sqrtCfOver2 = uTau/max(magUp[facei],smallNum);
		scalar Stk = 1.42*pow(max(KsPlus,smallNum),-0.45)*pow(max(Pr,smallNum),-0.8);
		alphaTurb[facei] = 1/Prt_;

		if(KsPlus > 5)
		{
			scalar heatAnalogy = 1/(Prt_+sqrtCfOver2/max(Stk,smallNum))*1/(Prt_*aFun);
			alphaTurb[facei] = heatAnalogy;
		}
	}

    //operator==(rhow*tnutw/Prt_);
	operator==(rhow*tnutw*alphaTurb);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void siboalphatWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
	os.writeKeyword("Ks") << Ks_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    siboalphatWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
