from pathlib import Path

import pyarrow as pa
from pyarrow import parquet

HERE = Path(__file__).parent


def test_export_without_register(con):
    tab = con.sql("""SELECT ST_GeomFromText('POINT (0 1)') as geom;""").to_arrow_table()
    assert tab.schema.field("geom").metadata is None


def test_basic_export(geoarrow_con):
    tab = geoarrow_con.sql(
        """SELECT ST_GeomFromText('POINT (0 1)') as geom;"""
    ).to_arrow_table()
    assert tab.schema.field("geom").metadata == {
        b"ARROW:extension:metadata": b"{}",
        b"ARROW:extension:name": b"geoarrow.wkb",
    }


def test_basic_import(geoarrow_con):
    tab = geoarrow_con.sql(
        """SELECT ST_GeomFromText('POINT (0 1)') as geom;"""
    ).to_arrow_table()
    assert tab.schema.field("geom").metadata == {
        b"ARROW:extension:metadata": b"{}",
        b"ARROW:extension:name": b"geoarrow.wkb",
    }


def test_roundtrip_segments(geoarrow_con):
    segments_file = HERE.parent / "data" / "segments.parquet"
    raw_table = parquet.read_table(segments_file)
    assert raw_table.schema.field("geometry").metadata is None

    field = pa.field(
        "geometry", pa.binary(), metadata={"ARROW:extension:name": "geoarrow.wkb"}
    )
    schema = pa.schema([field])
    geo_table = pa.table([raw_table["geometry"]], schema=schema)

    geoarrow_table_wkt = geoarrow_con.sql(
        """SELECT ST_AsText(geometry) as wkt FROM geo_table"""
    ).to_arrow_table()
    geoparquet_table_wkt = geoarrow_con.sql(
        f"""SELECT ST_AsText(geometry) as wkt FROM "{segments_file}";""",
    ).to_arrow_table()
    assert geoarrow_table_wkt == geoparquet_table_wkt

    # While we're here, check the roundtrip output
    assert geoarrow_con.sql("""SELECT * from geo_table""").to_arrow_table() == geo_table
