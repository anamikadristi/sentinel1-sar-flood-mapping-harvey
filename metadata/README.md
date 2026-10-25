# Metadata

FGDC Content Standard for Digital Geospatial Metadata (CSDGM, FGDC-STD-001-1998) records for the main geospatial products of this repository.

| Record | Data file | Product |
|---|---|---|
| [`flood_classes_20170830_20m_fgdc.xml`](flood_classes_20170830_20m_fgdc.xml) | `outputs/gis/flood_classes_20170830_20m.tif` | Sentinel-1 flood classes for 30 August 2017 (GeoTIFF, 20 m) |
| [`harvey_sar_flood_2017_fgdc.xml`](harvey_sar_flood_2017_fgdc.xml) | `outputs/gis/harvey_sar_flood_2017.gpkg` | Tract and ZCTA flood exposure and flood extent polygons (GeoPackage, three layers) |

## What each record contains

Each XML file has the full CSDGM sections: identification (citation, abstract, purpose, time period, status, bounding coordinates, theme and place keywords with thesauri, access and use constraints), data quality (attribute accuracy, logical consistency, completeness, positional accuracy, lineage with source citations and processing steps), spatial data organization, spatial reference, entity and attribute definitions, distribution and metadata reference. Extents, coordinate systems, cell sizes, value ranges, class counts and field lists were read from the data files themselves; methods and accuracy figures come from the scripts and the main README.

## Validation

* Every file is well formed XML.
* Every file validates against the official FGDC CSDGM XML Schema (`fgdc-std-001-1998.xsd` and its section modules, downloaded from https://www.fgdc.gov/schemas/metadata/).
* Every file passes the USGS metadata parser `mp` version 2.9.50 (built from the USGS source at https://geology.usgs.gov/tools/metadata/) with no errors.

The standard requires a mailing address and a voice telephone for each contact. Only the city and state are given; the postal code and telephone fields read "Not provided". Please use email for contact.

## Viewing

Open the XML in a browser or a metadata editor, or turn it into text or HTML with `mp`, for example `mp flood_classes_20170830_20m_fgdc.xml -t record.txt`.

Contact: Anamika Dristi, anamikadristi@gmail.com
