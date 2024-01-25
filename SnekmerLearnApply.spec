/*
A KBase module: SnekmerLearnApply
This sample module contains one small method that filters contigs.
*/

module SnekmerLearnApply {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    funcdef run_SnekmerLearnApply(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;




};
