/*
 * main.cpp
 *
 *  Created on: Jul 11, 2023
 *      Author: pemueller
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Core>
#include <vector>
#include <math.h>
#include <algorithm>
#include <boost/make_shared.hpp>
#include <cairo-pdf.h>
#include <cairo-svg.h>
#include <stdio.h>

#include "cellParameters.hpp"

using namespace std;

void readCellParameters( fstream& file, UnitCellPtr unitCell )
{
	file.seekg(0);
	string line, s;
	stringstream ss;

	while ( s != "CELLP" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	getline(file, line); // line now contains the cell parameters
	Vector3d lengths, angles;
	ss.clear();
	ss.str("");
	ss << line;

	Matrix3d lattice = Matrix3d::Zero();

	ss >> lattice(0,0);
	ss >> lattice(1,1);
	ss >> lattice(2,2);
	ss >> angles(0);
	ss >> angles(1);
	ss >> angles(2);

	for ( int i = 0; i < 3; ++i )
	{
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		double x = -1 * lattice(k,0) / lattice.row(k).norm();
		double y = -1 * lattice(k,1) / lattice.row(k).norm();
		double z = -1 * lattice(k,2) / lattice.row(k).norm();

		double cos_a = cos( M_PI/180 * ( angles(k) - 90 ) );
		double sin_a = sin( M_PI/180 * ( angles(k) - 90 ) );

		Matrix3d rotationMatrix;

		rotationMatrix(0,0) = cos_a + x * x * (1 - cos_a);
		rotationMatrix(0,1) = -1.0 * z * sin_a + x * y * (1 - cos_a);
		rotationMatrix(0,2) = y * sin_a + x * z * (1 - cos_a);
		rotationMatrix(1,0) = z * sin_a + x * y * (1 - cos_a);
		rotationMatrix(1,1) = cos_a + y * y * (1 - cos_a);
		rotationMatrix(1,2) = -1.0 * x * sin_a + y * z * (1 - cos_a);
		rotationMatrix(2,0) = -1.0 * y * sin_a + z * x * (1 - cos_a);
		rotationMatrix(2,1) = x * sin_a + z * y * (1 - cos_a);
		rotationMatrix(2,2) = cos_a + z * z * (1 - cos_a);

		lattice.row(j) = lattice.row(j) * rotationMatrix; // rotate the basis vectors
	}

	unitCell->setLattice(lattice);
}

void readBoundary( fstream& file, UnitCellPtr unitCell )
{
	file.seekg(0);
	string line, s;
	stringstream ss;

	while ( s != "BOUND" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}
	getline(file, line);
	ss.clear();
	ss.str("");
	ss << line;

	Vector3d from;
	Vector3d to;
	for ( int i = 0; i < 3; ++i )
	{
		ss >> from(i);
		ss >> to(i);
	}

	unitCell->setBoundary(make_pair(from,to));
}

void readAtoms( fstream& file, UnitCellPtr unitCell )
{
	file.seekg(0);
	string line, s;
	stringstream ss;
	const pair<Vector3d,Vector3d> boundary = unitCell->getBoundary();

	while ( s != "STRUC" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	while( true )
	{
		stringstream ss;
		string s;

		getline(file, line);

		ss.clear();
		ss.str("");
		ss << line;

		ss >> s;
		if ( s == "0" )
			break;

		string element;
		Vector3d coordinates;
		ss >> element;
		ss >> s >> s;
		for ( int i = 0; i < 3; ++i )
			ss >> coordinates(i);

		// apply boundary
		for ( int x = floor(boundary.first(0)); x <= ceil(boundary.second(0)); ++x )
		{
			for ( int y = floor(boundary.first(1)); y <= ceil(boundary.second(1)); ++y )
			{
				for ( int z = floor(boundary.first(2)); z <= ceil(boundary.second(2)); ++z )
				{
					Vector3d translation((double)x,(double)y,(double)z);
					Vector3d coordinatesNew  = coordinates + translation;

					if ( !( (coordinatesNew - boundary.first).minCoeff() < 0 || (coordinatesNew - boundary.second).maxCoeff() > 0) )
					{
						unitCell->addAtom( boost::make_shared<Atom>(element, coordinatesNew, unitCell->getLattice().transpose() * coordinatesNew) );
					}
				}
			}
		}

		getline(file, line);
	}
}

void readAtomAppearance( fstream& file, UnitCellPtr unitCell )
{
	file.seekg(0);
	string line, s, element;
	stringstream ss;

	while ( s != "ATOMT" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	while ( true )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;

		ss >> s;
		if ( s == "0" )
			break;

		ss >> element;
		double size;
		Vector3d rgb;
		ss >> size;
		for ( int i = 0; i < 3; ++i )
			ss >> rgb(i);

		for( size_t i = 0; i < unitCell->getAtoms().size(); ++i )
		{
			if ( unitCell->getAtoms()[i]->getElement() == element )
			{
				unitCell->getAtoms()[i]->setSize(size);
				unitCell->getAtoms()[i]->setRGB(rgb);
			}
		}
	}
}

void bondGeneratorUnitCell( UnitCellPtr unitCell, pair<string,string> types, pair<double,double> length, Vector3d rgb, double size )
{
	vector <AtomPtr> atoms = unitCell->getAtoms();

	for ( size_t mu = 0; mu < atoms.size(); ++mu )
	{
		AtomPtr atomMu = atoms[mu];

		if( types.first != atomMu->getElement() )
			continue;

		for ( size_t nu = mu + 1; nu < atoms.size(); ++nu )
		{
			AtomPtr atomNu = atoms[nu];

			if( types.second != atomNu->getElement() )
				continue;

			double d = (atomNu->getCoordinatesReal() - atomMu->getCoordinatesReal()).norm();
			if ( d < length.second && d > length.first  )
				unitCell->addBond(boost::make_shared<Bond>( make_pair(atomMu,atomNu), size, rgb ));
		}
	}
}

void bondGeneratorTranslation( UnitCellPtr unitCell, pair<string,string> types, pair<double,double> length, Vector3d rgb, double size )
{
	vector <AtomPtr> atoms = unitCell->getAtoms();

	Vector3d max;
	for ( size_t i = 0; i < 3; ++i )
		max(i) = ceil( length.second / unitCell->getLattice().row(i).norm() );


	for ( size_t mu = 0; mu < atoms.size(); ++mu )
	{
		AtomPtr atomMu = atoms[mu];

		if( !( types.first == atomMu->getElement() ))//|| types.second == atomMu->getElement() ))
			continue;

		for ( size_t nu = 0; nu < atoms.size(); ++nu )
		{
			AtomPtr atomNu = atoms[nu];

			if( !(( types.second == atomNu->getElement() )))//&& types.first == atomMu->getElement()) || ( types.second == atomMu->getElement() && types.first == atomNu->getElement() ) ))
				continue;

			// apply boundary to second atom
			for ( int x = -max(0); x <= max(0); ++x )
			{
				for ( int y = -max(1); y <= max(1); ++y )
				{
					for ( int z = -max(2); z <= max(2); ++z )
					{
						Vector3d translation((double)x,(double)y,(double)z);

						translation = translation + atomMu->getTranslation();

						AtomPtr atomNuT = boost::make_shared<Atom>( atomNu->getElement(), atomNu->getCoordinatesFrac() + translation, atomNu->getCoordinatesReal() + unitCell->getLattice().transpose() * translation, atomNu->getRGB() );

						double d = (atomNuT->getCoordinatesReal() - atomMu->getCoordinatesReal()).norm();

						if ( d < length.second && d > length.first )
						{
							if ( !unitCell->atomExists(atomNuT) )
							{
								//cout << "adding atom " << atomNuT->getElement() << " on " << atomNuT->getCoordinatesReal().transpose() << endl;
								unitCell->addAtom( atomNuT );
								unitCell->addBond(boost::make_shared<Bond>( make_pair(atomMu,atomNuT), size, rgb ));
							}
							else
							{
								AtomPtr atomNuT2 = unitCell->getAtom(atomNuT);
                                                                //cout << "atom " << atomNuT2->getElement() << " on " << atomNuT->getCoordinatesReal().transpose() << " already exists" << endl;

								unitCell->addBond(boost::make_shared<Bond>( make_pair(atomMu,atomNuT2), size, rgb ));
							}
						}
					}
				}
			}
		}
	}
}



void bondGenerator( UnitCellPtr unitCell, pair<string,string> types, pair<double,double> length, Vector3d rgb, double size, int boundary )
{
	if ( boundary == 0 ) // no translation at all
		bondGeneratorUnitCell(unitCell, types, length, rgb, size);

	if ( boundary == 1 ) // simple translation
	{
		//bondGeneratorUnitCell(unitCell, types, length, rgb, size);
		bondGeneratorTranslation(unitCell, types, length, rgb, size);
	}

	if ( boundary == 2 ) // recursive search
	{
		bondGeneratorUnitCell(unitCell, types, length, rgb, size);
		size_t oldSize = 0;
		do
		{
			cout << "# of atoms " << unitCell->getAtoms().size() << endl;

			oldSize = unitCell->getAtoms().size();
			bondGeneratorTranslation(unitCell, types, length, rgb, size);
		}

		while( oldSize != unitCell->getAtoms().size() );
		cout << "# of atoms after search " << unitCell->getAtoms().size() << endl;
	}
}

void convertPresetsToBonds(UnitCellPtr unitCell)
{
	for ( auto preset : unitCell->getBondPresets() )
	{
		if ( preset->getBoundary() != 0 )
			continue;

		bondGeneratorUnitCell( unitCell, preset->getTypes(), preset->getLength(), preset->getRGB(), preset->getSize() );
	}

	for ( auto preset : unitCell->getBondPresets() )
	{
		if ( preset->getBoundary() != 1 )
			continue;
		bondGeneratorTranslation(unitCell, preset->getTypes(), preset->getLength(), preset->getRGB(), preset->getSize());
	}

	int oldSize = 0;
	do
	{
		oldSize = unitCell->getAtoms().size();
		for ( auto preset : unitCell->getBondPresets() )
		{
			if ( preset->getBoundary() != 2 )
				continue;

			bondGeneratorTranslation(unitCell, preset->getTypes(), preset->getLength(), preset->getRGB(), preset->getSize());
			bondGeneratorTranslation(unitCell, make_pair(preset->getTypes().second, preset->getTypes().first), preset->getLength(), preset->getRGB(), preset->getSize());
		}
	}while( oldSize != unitCell->getAtoms().size() );
}

void readBonds( fstream& file, UnitCellPtr unitCell )
{
	file.seekg(0);

	// generate bonds
	string line, s;
	stringstream ss;

	while ( s != "SBOND" && file.good())
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}
	while( true )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;

		ss >> s;
		if ( s == "0" )
			break;

		pair<string,string> types;
		pair<double,double> length;
		double size;
		Vector3d rgb;
		int boundary;

		ss >> types.first;
		ss >> types.second;
		ss >> length.first;
		ss >> length.second;
		ss >> s;
		ss >> boundary;
		for ( int i = 0; i < 4; ++i )
			ss >> s;
		ss >> size;
		for ( int i = 0; i < 3; ++i )
			ss >> rgb(i);

		// make this into a template and do the search afterwards
		unitCell->addBondPreset( boost::make_shared<BondPreset>(types, length, rgb, size, boundary) );
	}

	convertPresetsToBonds(unitCell);
}

Vector3d getAtomFromFile( fstream& file, UnitCellPtr unitCell, int index )
{
	const auto prevPos = file.tellg();

	Vector3d result;

	string line, s;
	stringstream ss;

	file.seekg(0);

	while ( s != "STRUC" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	while ( true )
	{
		int i;
		getline(file,line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> i;

		if ( i == index )
		{
			ss >> s >> s >> s;

			for ( int j = 0; j < 3; ++j )
				ss >> result(j);

			result = unitCell->getLattice().transpose() * result;

			break;
		}

		getline(file,line);
	}

	file.seekg(prevPos);
	return result;
}

void readBondsAsVectors( fstream& file, UnitCellPtr unitCell )
{
	// Read rgb values and sizes of the vectors
	vector<pair<double, Vector3d>> props;
	file.seekg(0);

	string line, s;
	stringstream ss;

	while ( s != "VECTT" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	while (true)
	{
		double size;
		Vector3d rgb;
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
		if ( s == "0" )
			break;

		ss >> size;
		size *= 16;
		for ( size_t i = 0; i < 3; ++i )
			ss >> rgb(i);
		props.push_back(make_pair(size,rgb));
	}

	// read the vectors
	file.seekg(0);
	while ( s != "VECTR" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	while ( true )
	{
		int index;
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> index;
		if ( index == 0 )
			break;

		--index;
		Vector3d direction;
		for ( size_t i = 0; i < 3; ++i )
			ss >> direction(i);

		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		int atomIndex;
		ss >> atomIndex;

		Vector3d T;
		ss >> s;
		for ( size_t i = 0; i < 3; ++i )
			ss >> T(i);

		AtomPtr start = unitCell->getAtomAtPosition(getAtomFromFile(file, unitCell, atomIndex) + unitCell->getLattice().transpose() * T);
		AtomPtr end = unitCell->getAtomAtPosition( start->getCoordinatesReal() + direction );
		unitCell->addBond(boost::make_shared<Bond>(make_pair( start,end ), props[index].first, props[index].second));

/*		cout << "adding a bond from " << start->getElement() << start->getIndex() << " on " << start->getCoordinatesReal().transpose() <<
				" to " << end->getElement() << end->getIndex() << " on " << end->getCoordinatesReal().transpose() << endl;
		cout << "props are " << props[index].first << " " << props[index].second.transpose() << endl;
*/
		getline(file,line);
	}
}

void readOrientation( fstream& file, UnitCellPtr unitCell )
{
	file.seekg(0);

	string line, s;
	stringstream ss;

	while ( s != "SCENE" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	Matrix3d orientation;

	for ( int i = 0; i < 3; ++i )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;

		Vector3d line;

		for ( int j = 0; j < 3; ++j )
		{
			double d;
			ss >> d;
			line(j) = d;
		}
		orientation.row(i) = line;
	}

	unitCell->setOrientation(orientation);
}

void addUnitCell( UnitCellPtr unitCell )
{
	for ( double x = 0; x <= 1; ++x )
	{
		for ( double y = 0; y <= 1; ++y )
		{
			for ( double z = 0; z <= 1; ++z )
			{
				Vector3d latticeVector (x,y,z);

				unitCell->addAtom( boost::make_shared<Atom>("NULL", latticeVector, unitCell->getLattice().transpose() * latticeVector));
			}
		}
	}

	for ( auto atom1 : unitCell->getAtoms() )
	{
		for ( auto atom2 : unitCell->getAtoms() )
		{
			if ( (atom1->getCoordinatesFrac() - atom2->getCoordinatesFrac()).norm() < 1.00001 )
			{
				unitCell->addBond(boost::make_shared<Bond>( make_pair(atom1,atom2), 1 ));
			}
		}
	}

}

void applyOrientation( UnitCellPtr unitCell )
{
	const Matrix3d orientation = unitCell->getOrientation();
	for ( auto atom : unitCell->getAtoms() )
	{
		Vector3d coordinatesRealNew = orientation * atom->getCoordinatesReal();

		atom->setCoordinatesFrac( Vector3d::Zero() );
		atom->setCoordinatesReal( coordinatesRealNew);
	}
}

bool PointerCompare ( const boost::shared_ptr<Object> l, const boost::shared_ptr<Object> r )
{
/*	if ( l->getZ() == r->getZ() )
	{
		return (l->getType() != "circle");
	}*/

	return l->getZ() < r->getZ();
}

enum BondDrawMode
{
	BOND_CENTER = 0,      // bonds drawn between atom centers (original behavior)
	BOND_CLIP_SCREEN = 1, // bonds clipped to the projected (2D) atom circle
	BOND_CLIP_REAL = 2    // bonds clipped to the real (3D) atom sphere, then projected
};

// Atom radii are drawn on screen as 0.4 * scale * size (see exportFile). Since the
// projection is orthographic, a sphere of real radius 0.4 * size projects to exactly
// that screen circle, so the same constant doubles as the "real" sphere radius.
double atomRadius( AtomPtr atom )
{
	return 0.4 * atom->getSize();
}

// Shrinks rA/rB proportionally if they would overlap/invert the segment between a and b.
template <typename VectorT>
void clampTrimRadii( const VectorT& a, const VectorT& b, double& rA, double& rB )
{
	double len = (b - a).norm();
	if ( rA + rB > len && rA + rB > 1e-9 )
	{
		double k = len / (rA + rB);
		rA *= k;
		rB *= k;
	}
}

pair<Vector2d,Vector2d> bondEndpoints2D( AtomPtr atomA, AtomPtr atomB, int bondDrawMode )
{
	Vector2d a2 = atomA->getCoordinatesReal().head(2);
	Vector2d b2 = atomB->getCoordinatesReal().head(2);

	if ( bondDrawMode == BOND_CENTER )
		return make_pair(a2, b2);

	if ( bondDrawMode == BOND_CLIP_SCREEN )
	{
		double len = (b2 - a2).norm();
		if ( len < 1e-9 )
			return make_pair(a2, b2);

		double rA = atomRadius(atomA);
		double rB = atomRadius(atomB);
		clampTrimRadii(a2, b2, rA, rB);

		Vector2d unit = (b2 - a2) / len;
		return make_pair( a2 + rA * unit, b2 - rB * unit );
	}

	// BOND_CLIP_REAL: trim along the true 3D bond axis by the real sphere radius,
	// then project to 2D, so foreshortened bonds appear to dive into the sphere
	// rather than stopping exactly on the screen circle's rim.
	Vector3d a3 = atomA->getCoordinatesReal();
	Vector3d b3 = atomB->getCoordinatesReal();

	double len = (b3 - a3).norm();
	if ( len < 1e-9 )
		return make_pair(a2, b2);

	double rA = atomRadius(atomA);
	double rB = atomRadius(atomB);
	clampTrimRadii(a3, b3, rA, rB);

	Vector3d unit = (b3 - a3) / len;
	Vector3d start = a3 + rA * unit;
	Vector3d end = b3 - rB * unit;

	return make_pair( start.head(2), end.head(2) );
}

vector<ObjectPtr> generateObjects( UnitCellPtr unitCell, int bondDrawMode )
{
	vector<ObjectPtr> result;
	for ( auto atom : unitCell->getAtoms() )
	{
		result.push_back( boost::make_shared<Object>( "circle", atom->getCoordinatesReal().head(2), atom->getRGB(), atom->getCoordinatesReal()(2), atom->getSize() ));
	}

	for ( auto bond : unitCell->getBonds() )
	{
		pair<Vector2d,Vector2d> endpoints = bondEndpoints2D( bond->getFrom(), bond->getTo(), bondDrawMode );

		// In BOND_CLIP_REAL, cap the bond with a foreshortened ellipse representing
		// its own circular cross section (a true cylinder end cap), rather than a
		// flat cut. capForeshorten is the bond axis's z-component (how much it points
		// toward/away from the viewer): 0 = side-on (cap degenerates to a flat line),
		// 1 = pointing straight at the viewer (cap is a full circle).
		double capForeshorten = -1;
		if ( bondDrawMode == BOND_CLIP_REAL )
		{
			Vector3d axis = bond->getTo()->getCoordinatesReal() - bond->getFrom()->getCoordinatesReal();
			double axisLen = axis.norm();
			capForeshorten = ( axisLen > 1e-9 ) ? fabs(axis(2) / axisLen) : 1.0;
		}

		result.push_back( boost::make_shared<Object>( "line", endpoints.first, bond->getRGB(), bond->getCoordinatesReal()(2), bond->getSize(), endpoints.second, capForeshorten));
	}

	// sort objects by z value

	sort(result.begin(),result.end(),PointerCompare);

	return result;
}

double getXLength(vector<ObjectPtr> objects)
{
	double min = 0;
	double max = 0;

	for ( auto it : objects )
	{
		if ( it->getType() == "circle" )
		{
			if ( it->getPosition()(0) > max )
				max = it->getPosition()(0);
			if ( it->getPosition()(0) < min )
				min = it->getPosition()(0);
		}
		if ( it->getType() == "line" )
		{
			double xStart = it->getPosition()(0);
			double xEnd =  it->getDirection()(0);
			double maxX = std::max(xStart,xEnd);
			double minX = std::min(xStart,xEnd);

			if ( maxX > max )
				max = maxX;
			if ( minX < min )
				min = minX;
		}
	}

	for ( auto object : objects )
	{
		object->addX( -1 * min );
	}

	return max - min;
}

double getYLength(vector<ObjectPtr> objects)
{
	double min = 0;
	double max = 0;

	for ( auto it : objects )
	{
		if ( it->getType() == "circle" )
		{
			if ( it->getPosition()(1) + it->getSize() > max )
				max = it->getPosition()(1);
			if ( it->getPosition()(1)  - it->getSize() < min )
				min = it->getPosition()(1);
		}
		if ( it->getType() == "line" )
		{
			double yStart = it->getPosition()(1);
			double yEnd = it->getDirection()(1);
			double maxY = std::max(yStart,yEnd);
			double minY = std::min(yStart,yEnd);

			if ( maxY > max )
				max = maxY;
			if ( minY < min )
				min = minY;
		}
	}

	for ( auto object : objects )
	{
		object->addY( -1 * min );
	}

	return max - min;
}

void exportFile(vector<ObjectPtr> objects, string filename)
{
	string outputFilename = filename.substr(0,filename.size()-5) + "pdf";
	double scale = 20;	//
	double xMargin = 50;
	double yMargin = 50;


	double xLength = ceil(scale * getXLength(objects));
	double yLength = ceil(scale * getYLength(objects));

    cairo_surface_t *surface = cairo_pdf_surface_create(outputFilename.c_str(), xLength + 2*xMargin, yLength + 2 * yMargin );
    cairo_t *cr = cairo_create(surface);

    cout << "Exporting " << outputFilename << endl;

	for ( auto object : objects )
	{
		if ( object->getType() == "circle" )
		{
			cairo_set_line_width (cr, 1);
			Vector3d rgb = object->getRGB();
//			cout << object->getSize() << endl;

//			cout << "c" << xMargin + scale*object->getPosition()(0) << " " << yMargin + scale*object->getPosition()(1) << " at " << object->getZ() <<endl;

			cairo_arc(cr, xMargin + scale*object->getPosition()(0), yLength + yMargin - scale*object->getPosition()(1), 0.4*scale*object->getSize(), 0, 2 * M_PI);

			cairo_set_source_rgb (cr, rgb(0)/255.0, rgb(1)/255.0, rgb(2)/255.0);
			cairo_fill_preserve(cr);
			cairo_set_source_rgb (cr, 0, 0,0);
			cairo_stroke (cr);
		}

		if ( object->getType() == "line" )
		{
			cairo_set_line_width (cr, object->getSize());

/*			cout << "l from " << xMargin + scale*object->getPosition()(0) << " " << yMargin + scale*object->getPosition()(1) << " at " << object->getZ() << endl;
			cout << "l to " << xMargin + scale*object->getDirection()(0) << " " << yMargin + scale*object->getDirection()(1) << endl;
			cout << "l total " << (scale*object->getDirection() - scale*object->getPosition()).norm() << endl;
*/			Vector3d rgb = object->getRGB();
			cairo_set_source_rgb (cr, rgb(0)/255.0, rgb(1)/255.0,rgb(2)/255.0);

			double x0 = xMargin + scale*object->getPosition()(0);
			double y0 = yMargin + yLength - scale*object->getPosition()(1);
			double x1 = xMargin + scale*object->getDirection()(0);
			double y1 = yMargin + yLength - scale*object->getDirection()(1);

			cairo_move_to (cr, x0, y0);
			//cairo_rel_line_to(cr, scale*object->getDirection()(0), scale*object->getDirection()(1));
			cairo_line_to (cr, x1, y1);
//			cairo_close_path (cr);
			cairo_stroke (cr);

			// perspective-correct cylinder end caps (BOND_CLIP_REAL only)
			double capForeshorten = object->getCapForeshorten();
			if ( capForeshorten > 1e-3 )
			{
				double capRadius = 0.5 * object->getSize();
				double theta = atan2(y1 - y0, x1 - x0);
				double endsX[2] = { x0, x1 };
				double endsY[2] = { y0, y1 };

				for ( int end = 0; end < 2; ++end )
				{
					cairo_save(cr);
					cairo_translate(cr, endsX[end], endsY[end]);
					cairo_rotate(cr, theta);
					cairo_scale(cr, capForeshorten, 1.0);
					cairo_arc(cr, 0, 0, capRadius, 0, 2 * M_PI);
					cairo_fill(cr);
					cairo_restore(cr);
				}
			}
		}
	}

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

bool isStructureP1( fstream& file )
{
	file.seekg(0);

	string line, s;
	stringstream ss;
	int groupNr;

	while ( s != "GROUP" )
	{
		getline(file, line);
		ss.clear();
		ss.str("");
		ss << line;
		ss >> s;
	}

	getline(file, line);
	ss.clear();
	ss.str("");
	ss << line;
	ss >> groupNr;

	if( groupNr != 1 )
		return false;

	return true;
}

void processFile( string filename, int bondDrawMode )
{
	UnitCellPtr _unitCell = boost::make_shared<UnitCell>();

	// open .vesta file
	fstream file(filename);
	if( file.fail() )
	{
		cout << "File " << filename << " not found..." << endl;
		return;
	}

	if ( !isStructureP1(file) )
	{
		cout << "File " << filename << " seems to contain symmetrized atom positions." << endl << "However, I can only digest structures with P1 symmetry. Please adjust." << endl;
		return;
	}

	// read unit cell
	readCellParameters( file, _unitCell );
	addUnitCell(_unitCell);

	// read boundary
	readBoundary(file, _unitCell);
//	cout << "Boundary" << endl << _unitCell->getBoundary().first.transpose() << endl << _unitCell->getBoundary().second.transpose() << endl;

	// read atom positions + colors
	readAtoms( file, _unitCell );
	readAtomAppearance( file, _unitCell );

/*	cout << _unitCell->getAtoms().size() <<  " atoms" << endl;
	for ( auto atom : _unitCell->getAtoms() )
		cout << atom->getElement() << " " << atom->getCoordinatesFrac().transpose() << endl;
*/
	// read bonds
	readBondsAsVectors( file, _unitCell );
/*	cout << _unitCell->getBonds().size() << " after first step" << endl;
	cout << "they are" << endl;
	for ( auto bond : _unitCell->getBonds() )
		cout << bond->getFrom()->getElement() << bond->getFrom()->getIndex() << " (" << bond->getFrom()->getCoordinatesReal().transpose() << ")--"
		     << bond->getTo()->getElement() << bond->getTo()->getIndex() << " (" << bond->getTo()->getCoordinatesReal().transpose() << ") = " << (bond->getTo()->getCoordinatesReal().transpose() - bond->getFrom()->getCoordinatesReal().transpose()).norm() << endl;
*/

	readBonds( file, _unitCell );
//	cout << _unitCell->getBonds().size() << " after second step" << endl;
	readAtomAppearance( file, _unitCell );

	// read orientation
	readOrientation( file, _unitCell );
	applyOrientation(_unitCell);

/*	cout << _unitCell->getBonds().size() << " bonds" << endl;
	for ( auto bond : _unitCell->getBonds() )
		cout << bond->getFrom()->getElement() << bond->getFrom()->getIndex() << " (" << bond->getFrom()->getCoordinatesReal().transpose() << ")--"
		     << bond->getTo()->getElement() << bond->getTo()->getIndex() << " (" << bond->getTo()->getCoordinatesReal().transpose() << ") = " << (bond->getTo()->getCoordinatesReal().transpose() - bond->getFrom()->getCoordinatesReal().transpose()).norm() << endl;

	cout << _unitCell->getAtoms().size() <<  " atoms after bond detection" << endl;
	for ( auto atom : _unitCell->getAtoms() )
		cout << atom->getElement() << " " << atom->getCoordinatesReal().transpose() << endl;


	cout << _unitCell->getAtoms().size() <<  " atoms after orientation change" << endl;
	for ( auto atom : _unitCell->getAtoms() )
		cout << atom->getElement() << " " << atom->getCoordinatesReal().transpose() << endl;
*/
	file.close();
	// add objects projected on the x/y plane of the screen
	vector<ObjectPtr> objects = generateObjects(_unitCell, bondDrawMode);

	// export the final file
	exportFile(objects, filename);
}

void printUsage()
{
	cout << "Please provide a .vesta file for me to process..." << endl;
	cout << "Optional: --bond-mode <0|1|2>" << endl;
	cout << "  0 = bonds drawn between atom centers (default)" << endl;
	cout << "  1 = bonds clipped to the projected atom circle (good for side-on views)" << endl;
	cout << "  2 = bonds clipped to the real atom sphere, perspective-correct as in VESTA" << endl;
}

int main( int argc, char** argv )
{
	int bondDrawMode = BOND_CENTER;
	vector<string> files;

	for ( int i = 1; i < argc; ++i )
	{
		string arg = argv[i];
		if ( arg == "--bond-mode" && i + 1 < argc )
			bondDrawMode = stoi(argv[++i]);
		else
			files.push_back(arg);
	}

	if ( files.empty() )
	{
		printUsage();
		return 1;
	}

	for ( auto& filename : files )
	{
		cout << "Processing " << filename << endl;
		processFile(filename, bondDrawMode);
	}

	 return 0;
}


