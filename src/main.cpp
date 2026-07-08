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
#include <map>
#include <math.h>
#include <algorithm>
#include <boost/make_shared.hpp>
#include <cairo-pdf.h>
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

// True z at the two already-trimmed drawn endpoints, found by linear interpolation
// between the atoms' real z values. p0/p1 always lie on the segment from atomA's to
// atomB's 2D position (trimming only moves points inward along that same line), so
// their fractional position along that line gives the right interpolation parameter -
// exact for BOND_CLIP_REAL, and a consistent approximation for the other two modes.
pair<double,double> bondEndpointZs( AtomPtr atomA, AtomPtr atomB, Vector2d p0, Vector2d p1 )
{
	Vector2d a2 = atomA->getCoordinatesReal().head(2);
	Vector2d b2 = atomB->getCoordinatesReal().head(2);
	double zA = atomA->getCoordinatesReal()(2);
	double zB = atomB->getCoordinatesReal()(2);

	double lenSq = (b2 - a2).squaredNorm();
	if ( lenSq < 1e-18 )
		return make_pair(zA, zB);

	double t0 = (p0 - a2).dot(b2 - a2) / lenSq;
	double t1 = (p1 - a2).dot(b2 - a2) / lenSq;

	return make_pair( zA + t0 * (zB - zA), zA + t1 * (zB - zA) );
}

struct BondSegment
{
	Vector2d p0, p1;
};

// Finds the t (in [0,1]) along segment [p0,p1] where it crosses the boundary of circle
// (center c, radius r), localized to just one endpoint's neighborhood - used to confine a
// segment's own-atom occlusion-key boost (see generateObjects) to only the portion that's
// actually inside that atom's 2D disc, instead of inflating the whole segment's key.
// fromStart=true checks whether p0 starts inside the circle and returns where it exits
// moving forward (the larger root); fromStart=false checks p1 and returns where it exits
// moving backward (the smaller root). Returns -1 if that endpoint isn't inside the circle
// at all (no localized boost needed there).
double findDiscBoundaryT( Vector2d p0, Vector2d p1, Vector2d c, double r, bool fromStart )
{
	Vector2d probe = fromStart ? p0 : p1;
	if ( (probe - c).norm() >= r )
		return -1;

	Vector2d d = p1 - p0;
	Vector2d f = p0 - c;

	double a = d.dot(d);
	if ( a < 1e-18 )
		return -1;

	double b = 2 * f.dot(d);
	double cc = f.dot(f) - r * r;
	double discriminant = b * b - 4 * a * cc;
	if ( discriminant < 0 )
		return -1;

	double sqrtDisc = sqrt(discriminant);
	double t = fromStart ? (-b + sqrtDisc) / (2 * a) : (-b - sqrtDisc) / (2 * a);
	return min(1.0, max(0.0, t));
}

// Finds where the bond's true 3D path (p0,z0) to (p1,z1) passes through the inside of
// atom's sphere (true 3D center, true 3D radius r), and appends that t-interval to
// dropIntervals - it's genuinely occluded there, being literally embedded in the
// sphere's opaque volume, not just behind its flat center plane or its 2D silhouette.
// This is an exact 3D line-sphere intersection, not a 2D circle test: a bond can dip
// inside a sphere and re-emerge partway through its silhouette (whenever the far
// endpoint is enough nearer the viewer to overtake the sphere's near surface before
// reaching the silhouette's rim), and only a true 3D test captures that partial split -
// a single front/behind sample over the whole 2D-overlap interval cannot. A bond's own
// endpoint atom falls out of this for free: right at it, the 3D path point is exactly
// the atom's own center, deepest inside its own sphere, so the interval starts there
// and only ends where the path actually exits the sphere.
void computeOcclusionInterval( Vector2d p0, Vector2d p1, double z0, double z1,
                                Vector3d center, double r,
                                vector<pair<double,double>>& dropIntervals )
{
	Vector3d p0_3(p0(0), p0(1), z0);
	Vector3d p1_3(p1(0), p1(1), z1);

	Vector3d d = p1_3 - p0_3;
	Vector3d f = p0_3 - center;

	double a = d.dot(d);
	if ( a < 1e-18 )
		return;

	double b = 2 * f.dot(d);
	double cc = f.dot(f) - r * r;
	double discriminant = b * b - 4 * a * cc;

	if ( discriminant < 0 )
		return; // the bond's 3D path never enters the sphere

	double sqrtDisc = sqrt(discriminant);
	double tEnter = (-b - sqrtDisc) / (2 * a);
	double tExit = (-b + sqrtDisc) / (2 * a);

	double tStart = max(0.0, tEnter);
	double tEnd = min(1.0, tExit);

	if ( tEnd - tStart < 1e-9 )
		return; // overlap with the drawn segment is empty (or a single point)

	dropIntervals.push_back( make_pair(tStart, tEnd) );
}

// Clips a bond's drawn segment against every atom it overlaps in 2D, removing the
// portions that are genuinely behind some atom (by true interpolated depth, not by
// sort-key approximation), and returns the surviving visible piece(s).
vector<BondSegment> clipBondAgainstAtoms( Vector2d p0, Vector2d p1, double z0, double z1, const vector<AtomPtr>& atoms )
{
	vector<pair<double,double>> dropIntervals;

	for ( auto atom : atoms )
	{
		double r = atomRadius(atom);
		if ( r < 1e-9 )
			continue;

		computeOcclusionInterval(p0, p1, z0, z1, atom->getCoordinatesReal(), r, dropIntervals);
	}

	if ( dropIntervals.empty() )
		return vector<BondSegment>{ BondSegment{p0, p1} };

	sort(dropIntervals.begin(), dropIntervals.end());

	vector<pair<double,double>> merged;
	for ( auto& interval : dropIntervals )
	{
		if ( !merged.empty() && interval.first <= merged.back().second + 1e-9 )
			merged.back().second = max(merged.back().second, interval.second);
		else
			merged.push_back(interval);
	}

	vector<BondSegment> result;
	double cursor = 0.0;
	for ( auto& interval : merged )
	{
		if ( interval.first - cursor > 1e-4 )
			result.push_back( BondSegment{ p0 + cursor * (p1 - p0), p0 + interval.first * (p1 - p0) } );
		cursor = max(cursor, interval.second);
	}
	if ( 1.0 - cursor > 1e-4 )
		result.push_back( BondSegment{ p0 + cursor * (p1 - p0), p1 } );

	return result;
}

vector<ObjectPtr> generateObjects( UnitCellPtr unitCell, int bondDrawMode, bool overrideBondColor = false, Vector3d bondColor = Vector3d::Zero(), double bondWidthScale = 1.0 )
{
	vector<ObjectPtr> result;
	for ( auto atom : unitCell->getAtoms() )
	{
		result.push_back( boost::make_shared<Object>( "circle", atom->getCoordinatesReal().head(2), atom->getRGB(), atom->getCoordinatesReal()(2), atom->getSize() ));
	}

	int nextBondId = 0;
	for ( auto bond : unitCell->getBonds() )
	{
		AtomPtr atomA = bond->getFrom();
		AtomPtr atomB = bond->getTo();
		// Wireframe unit-cell edges are also Bonds internally, connecting placeholder
		// "NULL" atoms (see addUnitCell) - exclude them from the color override so it
		// only ever touches real atom-to-atom bonds, matching the wireframe's
		// long-standing out-of-scope status for any rendering changes here.
		bool isWireframe = ( atomA->getElement() == "NULL" || atomB->getElement() == "NULL" );
		Vector3d bondRGB  = ( overrideBondColor && !isWireframe ) ? bondColor : bond->getRGB();
		double bondSize   = isWireframe ? bond->getSize() : bond->getSize() * bondWidthScale;
		pair<Vector2d,Vector2d> endpoints = bondEndpoints2D( atomA, atomB, bondDrawMode );

		// In BOND_CLIP_REAL, cap the bond with a foreshortened ellipse representing
		// its own circular cross section (a true cylinder end cap), rather than a
		// flat cut. capForeshorten is the bond axis's z-component (how much it points
		// toward/away from the viewer): 0 = side-on (cap degenerates to a flat line),
		// 1 = pointing straight at the viewer (cap is a full circle).
		double capForeshorten = -1;
		if ( bondDrawMode == BOND_CLIP_REAL )
		{
			Vector3d axis = atomB->getCoordinatesReal() - atomA->getCoordinatesReal();
			double axisLen = axis.norm();
			capForeshorten = ( axisLen > 1e-9 ) ? fabs(axis(2) / axisLen) : 1.0;
		}

		// Clip the bond's body against its own two atoms only, by true interpolated
		// depth - not by sort-key approximation - removing the portion genuinely behind
		// its far endpoint (the self-occlusion bug). Deliberately NOT clipped against
		// other, unrelated atoms: that would interrupt bonds in places the old z-key
		// approach never did (third-party ordering is out of scope and was working
		// acceptably as-is). Caps only go on whichever surviving sub-segment end is the
		// outermost (i.e. genuinely touches the bond's real atom there); internal clip
		// boundaries never get one.
		pair<double,double> zs = bondEndpointZs( atomA, atomB, endpoints.first, endpoints.second );
		vector<AtomPtr> ownAtoms = { atomA, atomB };
		vector<BondSegment> segments = clipBondAgainstAtoms( endpoints.first, endpoints.second, zs.first, zs.second, ownAtoms );

		for ( size_t i = 0; i < segments.size(); ++i )
		{
			double capStart = ( i == 0 ) ? capForeshorten : -1;
			double capEnd = ( i == segments.size() - 1 ) ? capForeshorten : -1;
			double baseZ = bond->getCoordinatesReal()(2);

			// A surviving segment's clip boundary against its own atom sits exactly where
			// the bond's 3D path crosses that atom's sphere surface - which, away from the
			// silhouette's dead center, is generally *inside* the atom's projected 2D disc
			// (the sphere bulges toward the viewer there). That sliver between the clip
			// boundary and the disc's rim is genuinely meant to be visible, in front of the
			// atom's surface - but the atom's own circle is a flat opaque fill covering its
			// whole disc, so it must be painted *before* this segment there, or it wrongly
			// hides that correctly-visible sliver. The shared bond sort key (near atom's
			// true z, unboosted) doesn't guarantee that. This affects BOND_CENTER (which
			// starts/ends literally at the atom's center) and also BOND_CLIP_REAL (whose
			// trim point sits exactly on the true 3D sphere surface, but - unless the trim
			// direction happens to be exactly perpendicular to the viewing axis - that
			// surface point's 2D *projection* still generally falls inside the atom's 2D
			// silhouette). BOND_CLIP_SCREEN is the only mode that's actually exempt: it
			// trims by the full 2D screen radius, so its endpoint sits exactly on the 2D
			// disc's rim, never inside it.
			//
			// The boost must be localized to just the portion still inside that atom's 2D
			// disc - boosting the *whole* segment (e.g. when one segment touches both atoms
			// and the two needed boosts differ) lets the bigger boost leak onto the part
			// near the *other* atom, where it's not needed, making that part draw after
			// unrelated third-party atoms/bonds it has no business being in front of.
			bool needA = false, needB = false; // does this segment need a localized boost for that atom at all?
			double tA = 0, tB = 1; // boundary of the boosted zone, only meaningful if needA/needB
			bool crossed = false; // protection zones overlap (short bond vs big atoms): fall back to one unsplit piece
			if ( bondDrawMode == BOND_CENTER || bondDrawMode == BOND_CLIP_REAL )
			{
				// Only the *farther* atom's end can ever legitimately show the bond poking
				// out in front of it (the path's depth, right after trimming/leaving that
				// atom's center, only rises above that atom's own center depth when the
				// *other* atom is nearer - otherwise it only ever sinks further behind,
				// staying occluded for as long as it remains within that atom's 2D disc).
				// Near the *nearer* atom's end, no boost is ever needed: that atom's own
				// circle already naturally draws on top via the plain unboosted key.
				bool aIsFarther = atomA->getCoordinatesReal()(2) < atomB->getCoordinatesReal()(2);
				bool bIsFarther = atomB->getCoordinatesReal()(2) < atomA->getCoordinatesReal()(2);

				// Splitting off a piece is only worth it if boosting would actually raise
				// the key above the bond's own plain baseZ - otherwise max(baseZ, boost)
				// reduces to baseZ everywhere anyway, and splitting just fragments the
				// bond into multiple identically-keyed pieces for no behavioural gain
				// (purely cosmetic harm for editing in a vector tool).
				if ( i == 0 && aIsFarther && atomA->getCoordinatesReal()(2) + atomA->getSize() + 1e-6 > baseZ )
				{
					double t = findDiscBoundaryT( segments[i].p0, segments[i].p1, atomA->getCoordinatesReal().head(2), atomRadius(atomA), true );
					if ( t >= 0 ) { needA = true; tA = t; }
				}
				if ( i == segments.size() - 1 && bIsFarther && atomB->getCoordinatesReal()(2) + atomB->getSize() + 1e-6 > baseZ )
				{
					double t = findDiscBoundaryT( segments[i].p0, segments[i].p1, atomB->getCoordinatesReal().head(2), atomRadius(atomB), false );
					if ( t >= 0 ) { needB = true; tB = t; }
				}
				if ( needA && needB && tA >= tB )
					crossed = true;
			}

			Vector2d p0 = segments[i].p0;
			Vector2d p1 = segments[i].p1;
			Vector2d dir = p1 - p0;

			vector<pair<double,double>> pieces; // (tStart, tEnd)
			if ( !crossed )
			{
				double cursor = 0.0;
				if ( needA ) { pieces.push_back({0.0, tA}); cursor = tA; }
				if ( needB ) { if ( tB > cursor ) pieces.push_back({cursor, tB}); cursor = tB; }
				if ( 1.0 - cursor > 1e-9 ) pieces.push_back({cursor, 1.0});
			}
			if ( pieces.empty() )
				pieces.push_back({0.0, 1.0});

			for ( size_t j = 0; j < pieces.size(); ++j )
			{
				double pStart = pieces[j].first, pEnd = pieces[j].second;
				double z = baseZ;
				if ( crossed )
				{
					// no clean unboosted middle to preserve anyway - boost for whichever
					// own atom(s) this segment end touches, same as a simple shared key.
					if ( needA ) z = max( z, atomA->getCoordinatesReal()(2) + atomA->getSize() + 1e-6 );
					if ( needB ) z = max( z, atomB->getCoordinatesReal()(2) + atomB->getSize() + 1e-6 );
				}
				else
				{
					if ( needA && pStart < tA + 1e-9 )
						z = max( z, atomA->getCoordinatesReal()(2) + atomA->getSize() + 1e-6 );
					if ( needB && pEnd > tB - 1e-9 )
						z = max( z, atomB->getCoordinatesReal()(2) + atomB->getSize() + 1e-6 );
				}

				double pcs = ( j == 0 ) ? capStart : -1;
				double pce = ( j == pieces.size() - 1 ) ? capEnd : -1;

				result.push_back( boost::make_shared<Object>( "line", p0 + pStart * dir, bondRGB, z, bondSize, p0 + pEnd * dir, pcs, pce, nextBondId));
			}
		}
		++nextBondId;
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
			double r = 0.4 * it->getSize(); // drawn radius in world units
			if ( it->getPosition()(0) + r > max ) max = it->getPosition()(0) + r;
			if ( it->getPosition()(0) - r < min ) min = it->getPosition()(0) - r;
		}
		if ( it->getType() == "line" )
		{
			double xStart = it->getPosition()(0);
			double xEnd =  it->getDirection()(0);
			double maxX = std::max(xStart,xEnd);
			double minX = std::min(xStart,xEnd);

			if ( maxX > max ) max = maxX;
			if ( minX < min ) min = minX;
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
			double r = 0.4 * it->getSize(); // drawn radius in world units
			if ( it->getPosition()(1) + r > max ) max = it->getPosition()(1) + r;
			if ( it->getPosition()(1) - r < min ) min = it->getPosition()(1) - r;
		}
		if ( it->getType() == "line" )
		{
			double yStart = it->getPosition()(1);
			double yEnd = it->getDirection()(1);
			double maxY = std::max(yStart,yEnd);
			double minY = std::min(yStart,yEnd);

			if ( maxY > max ) max = maxY;
			if ( minY < min ) min = minY;
		}
	}

	for ( auto object : objects )
	{
		object->addY( -1 * min );
	}

	return max - min;
}

// Returns the directory component of a path, or "." if there is none.
string dirName( const string& path )
{
	size_t sep = path.rfind('/');
	return ( sep == string::npos ) ? "." : path.substr(0, sep);
}

// Strip directory components and the trailing .vesta extension from a file path,
// returning just the bare name used to construct output filenames.
string baseName( const string& path )
{
	size_t sep = path.rfind('/');
	string name = ( sep == string::npos ) ? path : path.substr(sep + 1);
	if ( name.size() >= 6 && name.substr(name.size()-6) == ".vesta" )
		name = name.substr(0, name.size()-6);
	return name;
}

void exportFile( vector<ObjectPtr> objects, string filename, string outputDir = ".", double scale = 20, double atomOutline = 1.0 )
{
	string outputFilename = outputDir + "/" + baseName(filename) + ".pdf";
	double xMargin = 50;
	double yMargin = 50;


	double xLength = ceil(scale * getXLength(objects));
	double yLength = ceil(scale * getYLength(objects));

    cairo_surface_t *surface = cairo_pdf_surface_create(outputFilename.c_str(), xLength + 2*xMargin, yLength + 2 * yMargin );
    cairo_t *cr = cairo_create(surface);

    cout << "Exporting " << outputFilename << endl;

	for ( size_t idx = 0; idx < objects.size(); ++idx )
	{
		ObjectPtr object = objects[idx];

		if ( object->getType() == "circle" )
		{
			cairo_set_line_width (cr, atomOutline * (scale / 20.0));
			Vector3d rgb = object->getRGB();

			cairo_arc(cr, xMargin + scale*object->getPosition()(0), yLength + yMargin - scale*object->getPosition()(1), 0.4*scale*object->getSize(), 0, 2 * M_PI);

			cairo_set_source_rgb (cr, rgb(0)/255.0, rgb(1)/255.0, rgb(2)/255.0);
			cairo_fill_preserve(cr);
			cairo_set_source_rgb (cr, 0, 0,0);
			cairo_stroke (cr);
			continue;
		}

		if ( object->getType() == "line" )
		{
			double lineWidth = object->getSize() * (scale / 20.0);
			cairo_set_line_width (cr, lineWidth);

			Vector3d rgb = object->getRGB();
			cairo_set_source_rgb (cr, rgb(0)/255.0, rgb(1)/255.0,rgb(2)/255.0);

			double x0 = xMargin + scale*object->getPosition()(0);
			double y0 = yMargin + yLength - scale*object->getPosition()(1);
			double x1 = xMargin + scale*object->getDirection()(0);
			double y1 = yMargin + yLength - scale*object->getDirection()(1);

			cairo_move_to (cr, x0, y0);
			cairo_line_to (cr, x1, y1);
			cairo_stroke (cr);

			// perspective-correct cylinder end caps (BOND_CLIP_REAL only). capStart/capEnd
			// are independent so an internal clip boundary (from occlusion splitting) never
			// gets a cap - only a sub-segment end that touches a genuine atom does.
			double capStart = object->getCapStart();
			double capEnd = object->getCapEnd();
			if ( capStart > 1e-3 || capEnd > 1e-3 )
			{
				double capRadius = 0.5 * object->getSize() * (scale / 20.0);
				double theta = atan2(y1 - y0, x1 - x0);

				if ( capStart > 1e-3 )
				{
					cairo_save(cr);
					cairo_translate(cr, x0, y0);
					cairo_rotate(cr, theta);
					cairo_scale(cr, capStart, 1.0);
					cairo_arc(cr, 0, 0, capRadius, 0, 2 * M_PI);
					cairo_fill(cr);
					cairo_restore(cr);
				}
				if ( capEnd > 1e-3 )
				{
					cairo_save(cr);
					cairo_translate(cr, x1, y1);
					cairo_rotate(cr, theta);
					cairo_scale(cr, capEnd, 1.0);
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

void exportSVG( vector<ObjectPtr> objects, string filename, string outputDir = ".", double scale = 20, double atomOutline = 1.0 )
{
	string outputFilename = outputDir + "/" + baseName(filename) + ".svg";
	double xMargin = 50;
	double yMargin = 50;

	double xLength = ceil(scale * getXLength(objects));
	double yLength = ceil(scale * getYLength(objects));

	double width  = xLength + 2 * xMargin;
	double height = yLength + 2 * yMargin;

	FILE* f = fopen(outputFilename.c_str(), "w");
	if ( !f )
	{
		cout << "Could not open " << outputFilename << " for writing." << endl;
		return;
	}

	cout << "Exporting " << outputFilename << endl;

	// Separate atom circles (bondId=-1) from bond segment/cap objects. Each bond's
	// objects are grouped into a <g> in the SVG so the bond can be recolored as a unit
	// in Illustrator. Bond bodies use a filled <path> (rectangle) rather than a stroked
	// <line> so that bond bodies and end caps share the same fill attribute - selecting
	// the group and changing fill updates the whole bond at once.
	//
	// Group sort key: max z within the group. This ensures a bond's perspective-boosted
	// sliver near its farther atom draws after that atom's circle (self-occlusion fix),
	// at the cost of very slightly approximate third-party ordering in rare edge cases -
	// an accepted trade-off consistent with the overall painter's-algorithm approach.
	vector<ObjectPtr> circles;
	std::map<int, vector<ObjectPtr>> bondGroupMap;

	for ( auto& obj : objects )
	{
		if ( obj->getBondId() < 0 )
			circles.push_back(obj);
		else
			bondGroupMap[obj->getBondId()].push_back(obj);
	}

	struct BondGroup { int id; double maxZ; };
	vector<BondGroup> bondGroups;
	for ( auto& kv : bondGroupMap )
	{
		double maxZ = -1e18;
		for ( auto& o : kv.second )
			maxZ = max(maxZ, o->getZ());
		bondGroups.push_back({ kv.first, maxZ });
	}
	sort(bondGroups.begin(), bondGroups.end(),
	     []( const BondGroup& a, const BondGroup& b ){ return a.maxZ < b.maxZ; });

	fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" "
	           "width=\"%.1f\" height=\"%.1f\" viewBox=\"0 0 %.1f %.1f\">\n",
	        width, height, width, height);

	// Merge-emit circles and bond groups in ascending z order
	size_t ci = 0, bg = 0;
	while ( ci < circles.size() || bg < bondGroups.size() )
	{
		double cz  = ( ci < circles.size()   ) ? circles[ci]->getZ()    : 1e18;
		double bgz = ( bg < bondGroups.size() ) ? bondGroups[bg].maxZ   : 1e18;

		if ( cz <= bgz )
		{
			// emit one atom circle
			auto& obj = circles[ci++];
			Vector3d rgb = obj->getRGB();
			double cx = xMargin + scale * obj->getPosition()(0);
			double cy = yLength + yMargin - scale * obj->getPosition()(1);
			double r  = 0.4 * scale * obj->getSize();
			fprintf(f, "  <circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" "
			           "fill=\"rgb(%d,%d,%d)\" stroke=\"black\" stroke-width=\"%.3f\"/>\n",
			        cx, cy, r,
			        (int)round(rgb(0)), (int)round(rgb(1)), (int)round(rgb(2)),
			        atomOutline * (scale / 20.0));
		}
		else
		{
			// emit one bond group
			int gid = bondGroups[bg++].id;
			auto& segs = bondGroupMap[gid];

			// All segments in a bond share the same color
			Vector3d rgb = segs[0]->getRGB();
			int ri = (int)round(rgb(0)), gi = (int)round(rgb(1)), bi = (int)round(rgb(2));

			// If any segment in the group has a perspective cap (mode 2), render bond
			// bodies as filled rectangles so bodies and caps share the fill attribute
			// and can be recolored together. Without caps (modes 0/1), use plain stroked
			// lines — simpler to manipulate in a vector editor.
			bool hasCaps = false;
			for ( auto& obj : segs )
				if ( obj->getCapStart() > 1e-3 || obj->getCapEnd() > 1e-3 )
					{ hasCaps = true; break; }

			fprintf(f, "  <g>\n");
			for ( auto& obj : segs )
			{
				double x0 = xMargin + scale * obj->getPosition()(0);
				double y0 = yLength + yMargin - scale * obj->getPosition()(1);
				double x1 = xMargin + scale * obj->getDirection()(0);
				double y1 = yLength + yMargin - scale * obj->getDirection()(1);

				if ( hasCaps )
				{
					// Filled rectangle so body and caps share the same fill attribute
					double ddx = x1 - x0, ddy = y1 - y0;
					double len = sqrt(ddx*ddx + ddy*ddy);
					if ( len > 1e-6 )
					{
						double nx = -ddy / len, ny = ddx / len;
						double h  = 0.5 * obj->getSize() * (scale / 20.0);
						double p0x = x0 + nx*h, p0y = y0 + ny*h;
						double p1x = x0 - nx*h, p1y = y0 - ny*h;
						double p2x = x1 - nx*h, p2y = y1 - ny*h;
						double p3x = x1 + nx*h, p3y = y1 + ny*h;
						fprintf(f, "    <path d=\"M %.3f,%.3f L %.3f,%.3f L %.3f,%.3f L %.3f,%.3f Z\" "
						           "fill=\"rgb(%d,%d,%d)\"/>\n",
						        p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y, ri, gi, bi);
					}

					// Perspective-correct end caps
					double capStart = obj->getCapStart();
					double capEnd   = obj->getCapEnd();
					if ( capStart > 1e-3 || capEnd > 1e-3 )
					{
						double capRadius = 0.5 * obj->getSize() * (scale / 20.0);
						double theta     = atan2(y1 - y0, x1 - x0);
						double thetaDeg  = theta * 180.0 / M_PI;

						if ( capStart > 1e-3 )
							fprintf(f, "    <ellipse cx=\"0\" cy=\"0\" rx=\"%.3f\" ry=\"%.3f\" "
							           "fill=\"rgb(%d,%d,%d)\" "
							           "transform=\"translate(%.3f,%.3f) rotate(%.3f)\"/>\n",
							        capRadius * capStart, capRadius, ri, gi, bi, x0, y0, thetaDeg);

						if ( capEnd > 1e-3 )
							fprintf(f, "    <ellipse cx=\"0\" cy=\"0\" rx=\"%.3f\" ry=\"%.3f\" "
							           "fill=\"rgb(%d,%d,%d)\" "
							           "transform=\"translate(%.3f,%.3f) rotate(%.3f)\"/>\n",
							        capRadius * capEnd, capRadius, ri, gi, bi, x1, y1, thetaDeg);
					}
				}
				else
				{
					// No caps: plain stroked line, simpler to work with in a vector editor
					fprintf(f, "    <line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" "
					           "stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%.3f\"/>\n",
					        x0, y0, x1, y1, ri, gi, bi, obj->getSize() * (scale / 20.0));
				}
			}
			fprintf(f, "  </g>\n");
		}
	}

	fprintf(f, "</svg>\n");
	fclose(f);
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

void processFile( string filename, int bondDrawMode, bool showCellEdges, bool outputSVG = false, string outputDir = "", double scale = 20, double bondWidthScale = 1.0, bool overrideBondColor = false, Vector3d bondColor = Vector3d::Zero(), double atomOutline = 1.0 )
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
	if ( showCellEdges )
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
	vector<ObjectPtr> objects = generateObjects(_unitCell, bondDrawMode, overrideBondColor, bondColor, bondWidthScale);

	// export the final file — default to alongside the input file
	string effectiveOutputDir = outputDir.empty() ? dirName(filename) : outputDir;
	if ( outputSVG )
		exportSVG(objects, filename, effectiveOutputDir, scale, atomOutline);
	else
		exportFile(objects, filename, effectiveOutputDir, scale, atomOutline);
}

void printUsage()
{
	cout << "Please provide a .vesta file for me to process..." << endl;
	cout << "Optional: --bond-mode <0|1|2>" << endl;
	cout << "  0 = bonds drawn between atom centers (default)" << endl;
	cout << "  1 = bonds clipped to the projected atom circle (good for side-on views)" << endl;
	cout << "  2 = bonds clipped to the real atom sphere, perspective-correct as in VESTA" << endl;
	cout << "Optional: --no-cell-edges" << endl;
	cout << "  Skip drawing the unit cell's wireframe box outline." << endl;
	cout << "Optional: --bond-color <name|R,G,B>" << endl;
	cout << "  Override every real bond's color (leaves the wireframe box untouched)." << endl;
	cout << "  Names: black, white, red, green, blue, gray. Or R,G,B with each 0-255." << endl;
	cout << "Optional: --pdf" << endl;
	cout << "  Output PDF instead of SVG (requires Cairo). SVG is the default." << endl;
	cout << "Optional: --svg" << endl;
	cout << "  Output SVG (default). Atoms import as single fill+stroke objects in" << endl;
	cout << "  Illustrator; bonds grouped for easy recoloring; no clipping masks." << endl;
	cout << "Optional: --output <directory>" << endl;
	cout << "  Directory to write the output file into (default: current directory)." << endl;
	cout << "Optional: --scale <number>" << endl;
	cout << "  Output scale factor in points (or SVG units) per Angstrom (default: 20)." << endl;
	cout << "  ~20 suits single-column figures; ~40-50 for larger or more detailed output." << endl;
	cout << "Optional: --bond-width <factor>" << endl;
	cout << "  Multiply all real bond widths by this factor (default: 1.5)." << endl;
	cout << "  Does not affect the unit cell wireframe." << endl;
	cout << "Optional: --atom-outline <factor>" << endl;
	cout << "  Atom outline stroke width as a multiple of scale/20 (default: 1.0)." << endl;
	cout << "  Use 0 to suppress outlines entirely." << endl;
}

// Parses --bond-color's argument: either a known name or an "R,G,B" triplet (0-255
// each). Returns false (and leaves color untouched) if the spec isn't recognized.
bool parseColor( const string& spec, Vector3d& color )
{
	if ( spec == "black" ) { color = Vector3d(0,0,0); return true; }
	if ( spec == "white" ) { color = Vector3d(255,255,255); return true; }
	if ( spec == "red" )   { color = Vector3d(255,0,0); return true; }
	if ( spec == "green" ) { color = Vector3d(0,255,0); return true; }
	if ( spec == "blue" )  { color = Vector3d(0,0,255); return true; }
	if ( spec == "gray" || spec == "grey" ) { color = Vector3d(128,128,128); return true; }

	stringstream ss(spec);
	string token;
	vector<double> components;
	while ( getline(ss, token, ',') )
	{
		try { components.push_back(stod(token)); }
		catch (...) { return false; }
	}
	if ( components.size() != 3 )
		return false;

	color = Vector3d(components[0], components[1], components[2]);
	return true;
}

int main( int argc, char** argv )
{
	int bondDrawMode = BOND_CENTER;
	bool showCellEdges = true;
	bool outputSVG = true;
	string outputDir = "";
	double scale = 20;
	double bondWidthScale = 1.5;
	double atomOutline = 1.0;
	bool overrideBondColor = false;
	Vector3d bondColor = Vector3d::Zero();
	vector<string> files;

	for ( int i = 1; i < argc; ++i )
	{
		string arg = argv[i];
		if ( arg == "--bond-mode" && i + 1 < argc )
			bondDrawMode = stoi(argv[++i]);
		else if ( arg == "--no-cell-edges" )
			showCellEdges = false;
		else if ( arg == "--svg" )
			outputSVG = true;
		else if ( arg == "--pdf" )
			outputSVG = false;
		else if ( arg == "--output" && i + 1 < argc )
			outputDir = argv[++i];
		else if ( arg == "--scale" && i + 1 < argc )
			scale = stod(argv[++i]);
		else if ( arg == "--bond-width" && i + 1 < argc )
			bondWidthScale = stod(argv[++i]);
		else if ( arg == "--atom-outline" && i + 1 < argc )
			atomOutline = stod(argv[++i]);
		else if ( arg == "--help" || arg == "-h" )
		{
			printUsage();
			return 0;
		}
		else if ( arg == "--bond-color" && i + 1 < argc )
		{
			string spec = argv[++i];
			if ( parseColor(spec, bondColor) )
				overrideBondColor = true;
			else
				cout << "Unrecognized --bond-color value \"" << spec << "\", ignoring." << endl;
		}
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
		processFile(filename, bondDrawMode, showCellEdges, outputSVG, outputDir, scale, bondWidthScale, overrideBondColor, bondColor, atomOutline);
	}

	 return 0;
}


