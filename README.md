A small Python library for doing conversions between Israeli grids (both Israel Transverse Mercator (ITM) aka רשת ישראל חדשה and Israel Cassini Soldner (ICS) aka רשת ישראל ישנה) and standard geographic coordinates.

Reimplemantes Joseph Gray's C++ library, read more [here](https://zvikabenhaim.appspot.com/software/ITM/).

### Example usage
```Python console
>>> import isr84
>>> isr84.wgs84_to_itm(32, 35)
(200131, 656329)
>>> isr84.wgs84_to_ics(32, 35)
(150145, 1156334)
>>> isr84.itm_to_wgs84(631556, 222286)
(31.7767479188685, 35.234383488170394)
```