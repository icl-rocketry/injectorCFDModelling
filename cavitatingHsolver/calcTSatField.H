{

	if (cryoFluid == "lh2")
	{
		const dimensionedScalar SatCoeff1("3.4812", dimensionSet(0,0,0,1,0,0,0), 3.4812);
		const dimensionedScalar SatCoeff2("19.82", dimensionSet(0,0,0,1,0,0,0), 19.82); 
		const dimensionedScalar plim("1", dimensionSet(1,-1,-2,0,0,0,0), 1.0); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		// from the code used in Fluent udf of MGDeGiorgi
		forAll(p, cellI)
		{
		if (p[cellI] <= plim.value())
		{ 
			TSat[cellI] = ( ( SatCoeff1.value() * Foam::log( plim.value()/plim.value() ) ) - SatCoeff2.value() );
		} else {
			TSat[cellI] = ( ( SatCoeff1.value() * Foam::log( p[cellI]/plim.value() ) ) - SatCoeff2.value() );
		}
		}
	
	} else if (cryoFluid == "ln2")
	{
		const dimensionedScalar SatCoeff1("8.1029", dimensionSet(0,0,0,1,0,0,0), 8.1029);
		const dimensionedScalar SatCoeff2("16.062", dimensionSet(0,0,0,1,0,0,0), 16.062); 
		const dimensionedScalar plim("1", dimensionSet(1,-1,-2,0,0,0,0), 1.0); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		// from the code used in Fluent udf of MGDeGiorgi
		forAll(p, cellI)
		{
		if (p[cellI] <= plim.value())
		{ 
			TSat[cellI] = ( ( SatCoeff1.value() * Foam::log( plim.value()/plim.value() ) ) - SatCoeff2.value() );
		} else {
			TSat[cellI] = ( ( SatCoeff1.value() * Foam::log( p[cellI]/plim.value() ) ) - SatCoeff2.value() );
		}
		}

	} else if (cryoFluid == "water")
	{	
		// from Affandi et al. (2012) - TSat(T) function
		Info << "Operating Fluid used: " << cryoFluid << endl;	
		const dimensionedScalar tCritical("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		//const dimensionedScalar tempMin("273.16", dimensionSet(0,0,0,1,0,0,0), 273.16);
		//const dimensionedScalar tempMax("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		const dimensionedScalar SatCoeff1("9.37817e-03", dimensionSet(0,0,0,0,0,0,0), 9.37817e-03);
		const dimensionedScalar SatCoeff2("4.98951e-04", dimensionSet(0,0,0,0,0,0,0), 4.98951e-04);
		const dimensionedScalar SatCoeff3("1.11049e-05", dimensionSet(0,0,0,0,0,0,0), 1.11049e-05);
		const dimensionedScalar SatCoeff4("3.34995e-07", dimensionSet(0,0,0,0,0,0,0), 3.34995e-07);
		const dimensionedScalar SatCoeff5("3.44102e-08", dimensionSet(0,0,0,0,0,0,0), 3.44102e-08);
		const dimensionedScalar ExpCoeff("-0.400", dimensionSet(0,0,0,0,0,0,0), -0.400);
		const dimensionedScalar pCritical("22064000", dimensionSet(1,-1,-2,0,0,0,0), 22064000); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		forAll(p, cellI)
		{
			TSat[cellI] = tCritical.value() * Foam::exp( ExpCoeff.value() * ( (SatCoeff1.value()) + ( SatCoeff2.value() * (p[cellI] / pCritical.value()) ) + ( SatCoeff3.value() * Foam::sqr(p[cellI] / pCritical.value()) ) + ( SatCoeff4.value() * Foam::pow( (p[cellI] / pCritical.value()),3.0) ) + ( SatCoeff5.value() * Foam::pow( (p[cellI] / pCritical.value()),4.0) ) ) );
		}

	} else {

		// By default water is considered -> Affandi et al. (2012) - TSat(T) function
		Info << "Operating Fluid used (default setting): water" << endl;	
		const dimensionedScalar tCritical("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		const dimensionedScalar SatCoeff1("9.37817e-03", dimensionSet(0,0,0,0,0,0,0), 9.37817e-03);
		const dimensionedScalar SatCoeff2("4.98951e-04", dimensionSet(0,0,0,0,0,0,0), 4.98951e-04);
		const dimensionedScalar SatCoeff3("1.11049e-05", dimensionSet(0,0,0,0,0,0,0), 1.11049e-05);
		const dimensionedScalar SatCoeff4("3.34995e-07", dimensionSet(0,0,0,0,0,0,0), 3.34995e-07);
		const dimensionedScalar SatCoeff5("3.44102e-08", dimensionSet(0,0,0,0,0,0,0), 3.44102e-08);
		const dimensionedScalar ExpCoeff("-0.400", dimensionSet(0,0,0,0,0,0,0), -0.400);
		const dimensionedScalar pCritical("22064000", dimensionSet(1,-1,-2,0,0,0,0), 22064000); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		forAll(p, cellI)
		{
			TSat[cellI] = tCritical.value() * Foam::exp( ExpCoeff.value() * ( (SatCoeff1.value()) + ( SatCoeff2.value() * (p[cellI] / pCritical.value()) ) + ( SatCoeff3.value() * Foam::sqr(p[cellI] / pCritical.value()) ) + ( SatCoeff4.value() * Foam::pow( (p[cellI] / pCritical.value()),3.0) ) + ( SatCoeff5.value() * Foam::pow( (p[cellI] / pCritical.value()),4.0) ) ) );
		}

	}

}
