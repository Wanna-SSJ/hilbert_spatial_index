#include <iostream>
#include <iomanip>
#include "hilbert.h"
#include "csv.h"
#include "btree.h"
#include <typeinfo>
using namespace std;

int main(int argv,char** argc)
{
    io::CSVReader<2> in("D:\\NYC\\yellow_tripdata_2016-01.csv");//green_tripdata_2016-01.csv"
    in.read_header(io::ignore_extra_column, "pickup_longitude","pickup_latitude");
    double pickup_lon,pickup_lat;

    stx::btree<uint64_t,uint32_t>* tree = new stx::btree<uint64_t,uint32_t>;
    int count = 0;
    while(in.read_row(pickup_lon,pickup_lat)){
    // do stuff with the data
        count ++;
        //cout << pickup_lon << "," << pickup_lat << endl;

        GeoPoint<GeoDegreeCoordinateSystem> point(pickup_lon,pickup_lat);

        GeoHasher<uint64_t> hasher;
        GeoHash<GeoDegreeCoordinateSystem,uint64_t> h= hasher.hash(point);
        tree->insert2(h.value(),count);
        //cout.precision(10);
        //cout << point.x() << endl;
        //cout << h.value() << endl;

        /*GeoPoint<GeoDegreeCoordinateSystem> point2 = hasher.unhash(h);
        cout << point2 << endl;
        break;*/
        if ( count % 100000 == 0)
            cout << tree->size() << endl;
    }

    cout << tree->size();
    GeoHasher<uint64_t> hasher;
    for (auto it = tree->begin(); it != tree->end(); it++) {
		//cout << typeid(uint64_t).name() << "," << it->second << endl;
		GeoHash<GeoDegreeCoordinateSystem,uint64_t> h(it->first);
        GeoPoint<GeoDegreeCoordinateSystem> point2 = hasher.unhash( h );
        cout.precision(10);
        //cout << point2 << endl;
	}

    return 0;

}
