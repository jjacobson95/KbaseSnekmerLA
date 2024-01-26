/*
A KBase module: SnekmerLearnApply
This sample module contains one small method that filters contigs.
This will have to be changed soon.
*/

module SnekmerLearnApply {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    funcdef run_SnekmerLearnApply(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;


};
