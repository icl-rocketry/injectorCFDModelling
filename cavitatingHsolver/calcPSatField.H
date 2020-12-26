{

	if (cryoFluid == "lh2")
	{
		//from Chung et al. 1976 - pSat(T) function
		Info << "Operating Fluid used: " << cryoFluid << endl;
		const dimensionedScalar tCritical("33.244", dimensionSet(0,0,0,1,0,0,0), 33.244);
		//const dimensionedScalar tempMin("13.95", dimensionSet(0,0,0,1,0,0,0), 13.95);
		//const dimensionedScalar tempMax("33.18", dimensionSet(0,0,0,1,0,0,0), 33.18);
		const dimensionedScalar SatCoeff1("2.65721", dimensionSet(0,0,0,0,0,0,0), 2.65721);
		const dimensionedScalar SatCoeff2("-2.74775", dimensionSet(0,0,0,0,0,0,0), -2.74775);
		const dimensionedScalar SatCoeff3("1.50557", dimensionSet(0,0,0,0,0,0,0), 1.50557);
		const dimensionedScalar SatCoeff4("0.10991", dimensionSet(0,0,0,0,0,0,0), 0.10991);
		const dimensionedScalar pSatMin("7175.41", dimensionSet(1,-1,-2,0,0,0,0), 7175.40); 
		const dimensionedScalar pSatMax("1309892.27", dimensionSet(1,-1,-2,0,0,0,0), 1309892.27); 
		const dimensionedScalar pCritical("1296656", dimensionSet(1,-1,-2,0,0,0,0), 1296656); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]


		forAll(T, cellI)
		{
			pSat[cellI] = pCritical.value() * Foam::exp( ( SatCoeff1.value() ) + ( (SatCoeff2.value() * tCritical.value()) / T[cellI] ) + ( SatCoeff3.value() * Foam::log( T[cellI]/tCritical.value() ) ) + ( SatCoeff4.value() * Foam::pow((T[cellI] / tCritical.value()),5.0) ) );
		}	

//pSat assumed constant
/*	
	forAll(T, cellI)
	{
		pSat[cellI] = 103444;
	}
*/

	} else if (cryoFluid == "ln2")
	{
		// from Span et al. 2000 - pSat(T) function
		Info << "Operating Fluid used: " << cryoFluid << endl;
		const dimensionedScalar tCritical("126.192", dimensionSet(0,0,0,1,0,0,0), 126.192);
		//const dimensionedScalar tempMin("63.151", dimensionSet(0,0,0,1,0,0,0), 63.151);
		//const dimensionedScalar tempMax("1000", dimensionSet(0,0,0,1,0,0,0), 1000);
		const dimensionedScalar SatCoeff1("-6.12445284", dimensionSet(0,0,0,0,0,0,0), -6.12445284);
		const dimensionedScalar SatCoeff2("1.26327220", dimensionSet(0,0,0,0,0,0,0), 1.26327220);
		const dimensionedScalar SatCoeff3("-0.765910082", dimensionSet(0,0,0,0,0,0,0), -0.765910082);
		const dimensionedScalar SatCoeff4("-1.77570564", dimensionSet(0,0,0,0,0,0,0), -1.77570564);
		const dimensionedScalar pCritical("3395800", dimensionSet(1,-1,-2,0,0,0,0), 3395800); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		forAll(T, cellI)
		{
			pSat[cellI] = pCritical.value() * Foam::exp( (tCritical.value()/T[cellI]) * ( ( SatCoeff1.value() * ( 1.0 - (T[cellI] / tCritical.value()) ) ) + ( SatCoeff2.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),1.5) ) + ( SatCoeff3.value()* Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),2.5) ) + ( SatCoeff4.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),5.0) ) ) );
		}	

	} else if (cryoFluid == "water")
	{
		//from Wagner and Pruss (1992) - pSat(T) function
		Info << "Operating Fluid used: " << cryoFluid << endl;	
		const dimensionedScalar tCritical("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		//const dimensionedScalar tempMin("273.16", dimensionSet(0,0,0,1,0,0,0), 273.16);
		//const dimensionedScalar tempMax("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		const dimensionedScalar SatCoeff1("-7.85951783", dimensionSet(0,0,0,0,0,0,0), -7.85951783);
		const dimensionedScalar SatCoeff2("1.84408259", dimensionSet(0,0,0,0,0,0,0), 1.84408259);
		const dimensionedScalar SatCoeff3("-11.7866497", dimensionSet(0,0,0,0,0,0,0), -11.7866497);
		const dimensionedScalar SatCoeff4("22.6807411", dimensionSet(0,0,0,0,0,0,0), 22.6807411);
		const dimensionedScalar SatCoeff5("-15.9618719", dimensionSet(0,0,0,0,0,0,0), -15.9618719);
		const dimensionedScalar SatCoeff6("1.80122502", dimensionSet(0,0,0,0,0,0,0), 1.80122502);
		const dimensionedScalar pCritical("22064000", dimensionSet(1,-1,-2,0,0,0,0), 22064000); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		forAll(T, cellI)
		{
			pSat[cellI] = pCritical.value() * Foam::exp( (tCritical.value()/T[cellI]) * ( ( SatCoeff1.value() * ( 1.0 - (T[cellI] / tCritical.value()) ) ) + ( SatCoeff2.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),1.5) ) + ( SatCoeff3.value()* Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),3.0) ) + ( SatCoeff4.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),3.5) ) + ( SatCoeff5.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),4.0) ) + ( SatCoeff6.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),7.5) ) ) );
		}

	} else {
		// by default water is considered as operating fluid -> Wagner and Pruss (1992) - pSat(T) function
		Info << "Operating Fluid used (default setting): water" << cryoFluid << endl;	
		const dimensionedScalar tCritical("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		//const dimensionedScalar tempMin("273.16", dimensionSet(0,0,0,1,0,0,0), 273.16);
		//const dimensionedScalar tempMax("647.096", dimensionSet(0,0,0,1,0,0,0), 647.096);
		const dimensionedScalar SatCoeff1("-7.85951783", dimensionSet(0,0,0,0,0,0,0), -7.85951783);
		const dimensionedScalar SatCoeff2("1.84408259", dimensionSet(0,0,0,0,0,0,0), 1.84408259);
		const dimensionedScalar SatCoeff3("-11.7866497", dimensionSet(0,0,0,0,0,0,0), -11.7866497);
		const dimensionedScalar SatCoeff4("22.6807411", dimensionSet(0,0,0,0,0,0,0), 22.6807411);
		const dimensionedScalar SatCoeff5("-15.9618719", dimensionSet(0,0,0,0,0,0,0), -15.9618719);
		const dimensionedScalar SatCoeff6("1.80122502", dimensionSet(0,0,0,0,0,0,0), 1.80122502);
		const dimensionedScalar pCritical("22064000", dimensionSet(1,-1,-2,0,0,0,0), 22064000); 

		// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*Sˆ2)]=[Pa]

		forAll(T, cellI)
		{
			pSat[cellI] = pCritical.value() * Foam::exp( (tCritical.value()/T[cellI]) * ( ( SatCoeff1.value() * ( 1.0 - (T[cellI] / tCritical.value()) ) ) + ( SatCoeff2.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),1.5) ) + ( SatCoeff3.value()* Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),3.0) ) + ( SatCoeff4.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),3.5) ) + ( SatCoeff5.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),4.0) ) + ( SatCoeff6.value() * Foam::pow( ( 1.0 - (T[cellI] / tCritical.value()) ),7.5) ) ) );
		}

	}
}
