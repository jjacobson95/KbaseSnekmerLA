# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except ImportError:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport


class cb_annotation_ontology_api(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login',
            service_ver='release',
            async_job_check_time_ms=100, async_job_check_time_scale_percent=150, 
            async_job_check_max_time_ms=300000):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = service_ver
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc,
            async_job_check_time_ms=async_job_check_time_ms,
            async_job_check_time_scale_percent=async_job_check_time_scale_percent,
            async_job_check_max_time_ms=async_job_check_max_time_ms)

    def get_annotation_ontology_events(self, params, context=None):
        """
        Retrieves annotation ontology events in a standardized form cleaning up inconsistencies in underlying data
        :param params: instance of type "GetAnnotationOntologyEventsParams"
           -> structure: parameter "input_ref" of String, parameter
           "input_workspace" of String, parameter "query_events" of list of
           String, parameter "query_genes" of list of String, parameter
           "standardize_modelseed_ids" of Long
        :returns: instance of type "GetAnnotationOntologyEventsOutput" ->
           structure: parameter "events" of list of type
           "AnnotationOntologyEvent" -> structure: parameter "event_id" of
           String, parameter "description" of String, parameter "ontology_id"
           of String, parameter "method" of String, parameter
           "method_version" of String, parameter "timestamp" of String,
           parameter "feature_types" of mapping from String to String,
           parameter "ontology_terms" of mapping from String to list of type
           "AnnotationOntologyTerm" -> structure: parameter "term" of String,
           parameter "modelseed_ids" of list of String, parameter
           "evidence_only" of Long, parameter "evidence" of type "Evidence"
           -> structure: parameter "reference" of tuple of size 2: parameter
           "entity_type" of String, parameter "ref_entity" of String,
           parameter "scores" of mapping from String to Double
        """
        return self._client.run_job('cb_annotation_ontology_api.get_annotation_ontology_events',
                                    [params], self._service_ver, context)

    def add_annotation_ontology_events(self, params, context=None):
        """
        Adds a new annotation ontology event to a genome or AMA
        :param params: instance of type "AddAnnotationOntologyEventsParams"
           -> structure: parameter "input_ref" of String, parameter
           "input_workspace" of String, parameter "output_name" of String,
           parameter "output_workspace" of String, parameter "clear_existing"
           of Long, parameter "overwrite_matching" of Long, parameter
           "events" of list of type "AnnotationOntologyEvent" -> structure:
           parameter "event_id" of String, parameter "description" of String,
           parameter "ontology_id" of String, parameter "method" of String,
           parameter "method_version" of String, parameter "timestamp" of
           String, parameter "feature_types" of mapping from String to
           String, parameter "ontology_terms" of mapping from String to list
           of type "AnnotationOntologyTerm" -> structure: parameter "term" of
           String, parameter "modelseed_ids" of list of String, parameter
           "evidence_only" of Long, parameter "evidence" of type "Evidence"
           -> structure: parameter "reference" of tuple of size 2: parameter
           "entity_type" of String, parameter "ref_entity" of String,
           parameter "scores" of mapping from String to Double
        :returns: instance of type "AddAnnotationOntologyEventsOutput" ->
           structure: parameter "output_ref" of String
        """
        return self._client.run_job('cb_annotation_ontology_api.add_annotation_ontology_events',
                                    [params], self._service_ver, context)

    def status(self, context=None):
        return self._client.run_job('cb_annotation_ontology_api.status',
                                    [], self._service_ver, context)
