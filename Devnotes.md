# Design decisions to consider

* Looping through multiple samples should be handled in the command line/SGE script (i.e. outside Bamboozle)
  * Aside from `barcode`, everything should work on a single sample
* Output file naming should reflect the command used, and then be compatible with utility scripts
* Implement more unit testing
  * Unit tests should be saved to `tests` directory
  * Example data for unit tests should be saved in `data/example_data`
