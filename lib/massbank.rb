#!/usr/bin/env ruby

require 'rubygems'
require 'rdf'
require 'rdf/ntriples'
require 'rdf/turtle'
require 'rdf/vocab'
include RDF

module MassBank

  class Record

    def initialize(txt)
      @accession
      @record = {}
      @pk_annot_head_1 = nil
      @pk_annot_head_2 = nil
      parse(txt)
    end
    attr_accessor :accession, :record, :mbo, :pk_annot_head_1, :pk_annot_head_2

    private

    def parse(txt)
      pk_annotation_flag = 0
      pk_peak_flag = 0

      txt.each_line do |line|
        line.chomp!
        case line
        when /PK\$ANNOTATION:\s+(.+)$/
          pk_annotation_flag = 1
          peak_annotation_header = $1
          if /^(.+)\s+\{(.+)\}/ =~ peak_annotation_header ||
             /^(.+)\s+\((.+)\)/ =~ peak_annotation_header
            annotation_line1 = $1
            annotation_line2 = $2
            @pk_annot_head_1 = annotation_line1.strip.chomp.split(/\s+/)
            @pk_annot_head_2 = annotation_line2.strip.chomp.split(/\s+/)
          else
            @pk_annot_head_1 = peak_annotation_header.strip.chomp.split(/\s+/)
          end
          @record["PK$ANNOTATION"] = []
        when /PK\$NUM_PEAK: (.+)/
          pk_annotation_flag = 0
          @record["PK$NUM_PEAK"] = [$1]
        when /PK\$PEAK/
          pk_peak_flag = 1
          @record["PK$PEAK"] = []
        when /\/\//
          pk_peak_flag = 0
        when /^\s\s(\S.+)$/
          pk_info = $1
          if pk_annotation_flag == 1
            @record["PK$ANNOTATION"] << Hash[@pk_annot_head_1.zip(pk_info.split(/\s+/))]
          elsif pk_peak_flag == 1
            @record["PK$PEAK"] << pk_info
          end
        when /^\s\s\s\s(\S.+)$/
          pk_info = $1
          if /^(\[lyso \S+)\s+(.+)$/ =~ pk_info
            pk_info_2 = [$1] + $2.split(/\s+/)
            @record["PK$ANNOTATION"][-1].update(Hash[@pk_annot_head_2.zip(pk_info_2)])
          else
            @record["PK$ANNOTATION"][-1].update(Hash[@pk_annot_head_2.zip(pk_info.split(/\s+/))])
          end
        when /^(\w.+?): (.+)$/
          @record.key?($1) ? @record[$1] << $2 : @record[$1] = [$2]
        when /^(\w.+):$/
        else
        end
      end
    end #end of parse
  end


  class RDFFactory

    def initialize(record, f_out)
      @record = record
      @accession = ""
      @section = ""
      @subtag = ""
      @f_out = f_out
      @statements = []
      @primary_subject
      @ch_bnode = nil
      # blank node for the Information of Chemical Compound Analyzed
      @sp_bnode = nil
      # blank node for the Information of Biological Sample
      @ac_bnode = nil
      # blank node for the Analytical Methods and Conditions
      @ms_bnode = nil
      # blank node for the Description of mass spectral data
      @obo = RDF::Vocabulary.new("http://purl.obolibrary.org/obo/")
      @sio = RDF::Vocabulary.new("http://semanticscience.org/resource/")
      @mbo = RDF::Vocabulary.new("http://www.massbank.jp/ontology/")
      @bibo = RDF::Vocabulary.new("http://purl.org/ontology/bibo/")
      @prefixes = {
        ns: "http://rdf.freebase.com/ns/",
        owl: "http://www.w3.org/2002/07/owl#>",
        rdfs: "http://www.w3.org/2000/01/rdf-schema#",
        rdf: "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        xsd: "http://www.w3.org/2001/XMLSchema#",
        dc: "http://purl.org/dc/terms/",
        skos: "http://www.w3.org/2004/02/skos/core#",
        obo: "http://purl.obolibrary.org/obo/",
        sio: "http://semanticscience.org/resource/",
        cas: "http://identifiers.org/cas/",
        chebi: "http://identifiers.org/chebi/",
        chempdb: "http://identifiers.org/pdb-ccd/",
        chemspider: "http://identifiers.org/chemspider/",
        inchikey: "http://identifiers.org/inchikey/",
        kegg: "http://identifiers.org/kegg.compound/",
        knapsack: "http://identifiers.org/knapsack/",
        nikkaji: "http://stirdf.jst.go.jp/cde/nikkaji/",
        lipidmaps: "http://identifiers.org/lipidmaps/",
        lipidbank: "http://lipidbank.jp/cgi-bin/detail.cgi?id=",
        pubchem: "http://identifiers.org/pubchem.substance/",
        pubchem_compound: "http://identifiers.org/pubchem.compound/",
        hmdb: "http://www.hmdb.ca/metabolites/",
        comptox: "http://comptox.epa.gov/dashboard/chemical/details/",
        taxonomy: "http://identifiers.org/taxonomy/",
        pubmed: "http://identifiers.org/pubmed/",
        bibo: "http://purl.org/ontology/bibo/",
        massbank: "http://www.massbank.jp/record/",
        mbo: "http://www.massbank.jp/ontology/"
      }
      @ion_type = ["[M]+",         "[M]+*",      "[M+H]+",         "[2M+H]+",
                   "[M+Na]+",      "[M-H+Na]+",  "[2M+Na]+",       "[M+2Na-H]+",
                   "[(M+NH3)+H]+", "[M+H-H2O]+", "[M+H-C6H10O4]+", "[M+H-C6H10O5]+",
                   "[M]-",         "[M-H]-",     "[M-2H]-",        "[M-2H+H2O]-",
                   "[M-H+OH]-",    "[2M-H]-",    "[M+HCOO-]-",     "[(M+CH3COOH)-H]-",
                   "[2M-H-CO2]-",  "[2M-H-C6H10O5]-"]
      @ion_type2id = Hash[@ion_type.zip((1..@ion_type.size).map{|i| "IT" + sprintf("%03d", i)})]
    end

    def rdfize
      @record.to_a.each do |key, value|
        v = value[0]
        @section = key
        case key
        when "ACCESSION"
          @statements << rdf_primary_resource(v)
          @statements << rdf_identifier(v)
        when "RECORD_TITLE"
          @statements << rdf_title(v)
        when "DATE"
          @statements = @statements + rdf_date(v)
        when "AUTHORS"
          @statements << rdf_authors(v)
        when "LICENSE"
          @statements << rdf_license(v)
        when "COPYRIGHT"
          @statements << rdf_copyright(v)
        when "PUBLICATION"
          @statements = @statements + rdf_publication(v)
        when "COMMENT"
          @statements << rdf_comment(v)
        when "CH$NAME"
          @statements = @statements + rdf_ch_name(value)
        when "CH$COMPOUND_CLASS"
          @statements << rdf_ch_compound_class(v)
        when "CH$FORMULA"
          @statements << rdf_ch_formula(v)
        when "CH$EXACT_MASS"
          @statements = @statements + rdf_ch_exact_mass(v)
        when "CH$SMILES"
          @statements << rdf_ch_smiles(v)
        when "CH$IUPAC"
          @statements << rdf_ch_iupac(v)
        when "CH$LINK"
          @statements = @statements + rdf_ch_link(value)
        when "SP$SCIENTIFIC_NAME"
          @statements = @statements + rdf_sp_scientific_name(v)
        when "SP$NAME"
          @statements << rdf_sp_name(v)
        when "SP$LINEAGE"
          @statements << rdf_sp_lineage(v)
        when "SP$LINK"
          @statements << rdf_sp_link(v)
        when "SP$SAMPLE"
          @statements << rdf_sp_sample(v)
        when "AC$INSTRUMENT"
          @statements = @statements + rdf_ac_instrument(v)
        when "AC$INSTRUMENT_TYPE"
          @statements << rdf_ac_instrument_type(v)
        when "AC$MASS_SPECTROMETRY"
          @statements = @statements + rdf_ac_mass_spectrometry(value)
        when "AC$CHROMATOGRAPHY"
          @statements = @statements + rdf_ac_chromatography(value)
        when "MS$FOCUSED_ION"
          @statements = @statements + rdf_ms_focused_ion(value)
        when "MS$DATA_PROCESSING"
          @statements << rdf_ms_data_processing(v)
        when "PK$SPLASH"
          @statements << rdf_pk_splash(v)
        when "PK$ANNOTATION"
          @statements = @statements + rdf_pk_annotation(value)
        when "PK$NUM_PEAK"
          @statements << rdf_pk_num_peak(v)
        when "PK$PEAK"
          @statements = @statements + rdf_pk_peak(value)
        else
        end
      end
      RDF::Turtle::Writer.open(@f_out, prefixes: @prefixes) do |writer|
        @statements.each do |statement|
          writer << statement
        end
      end
    end

    def rdf_primary_resource(v)
      @primary_subject = RDF::URI.new(@prefixes[:massbank] + v)
      @accession = v
      statement(@primary_subject, RDF.type, @mbo[:MassSpectrum])
    end

    def rdf_identifier(v)
      statement(@primary_subject, RDF::Vocab::DC.identifier, RDF::Literal.new(v))
    end

    def rdf_title(v)
      statement(@primary_subject,
                RDF::Vocab::SKOS.definition,
                RDF::Literal.new(v, :language => :en))
    end

    def rdf_date(v)
      statements = []
      date_current = ""
      date_created = ""
      date_modified = ""

      if /([\d\.]+)\s+\(Created\s+([\d\.]+),\s+modified\s+([\d\.]+)\)/ =~ v
        date_current = $1
        date_created = $2
        date_modified = $3

        statements << statement(@primary_subject,
                                RDF::Vocab::DC.date,
                                RDF::Literal.new(date_current.gsub(".", "-"),
                                                 :datatype => RDF::XSD.date))
        statements << statement(@primary_subject,
                                RDF::Vocab::DC.created,
                                RDF::Literal.new(date_current.gsub(".", "-"),
                                                    :datatype => RDF::XSD.date))
        statements << statement(@primary_subject,
                                RDF::Vocab::DC.modified,
                                RDF::Literal.new(date_current.gsub(".", "-"),
                                                 :datatype => RDF::XSD.date))
      elsif /([\d\.]+)/ =~ v
        date_current = $1
        statements << statement(@primary_subject,
                                RDF::Vocab::DC.date,
                                RDF::Literal.new(date_current.gsub(".", "-"),
                                                 :datatype => RDF::XSD.date))
      else
        put_error_message("Unknown DATE format .")
      end
      statements
    end

    def rdf_authors(v)
      statement(@primary_subject, RDF::Vocab::DC.creator, RDF::Literal.new(v))
    end

    def rdf_license(v)
      statement(@primary_subject, RDF::Vocab::DC.license, RDF::Literal.new(v))
    end

    def rdf_copyright(v)
      statement(@primary_subject, RDF::Vocab::DC.rights, RDF::Literal.new(v))
    end

    def rdf_publication(v)
      statements = []
      bnode = RDF::Node.new()
      statements << statement(@primary_subject, RDF::Vocab::DC.references, bnode)
      statements << statement(bnode, RDF.type, @bibo[:Article])
      statements << statement(bnode, RDF::RDFS.label, v)
      if /PMID:\s*(\d+)/ =~ v
        pmid = $1
        statements << statement(bnode, RDFS.seeAlso, RDF::URI.new(@prefixes[:pubmed] + pmid))
      end
      statements
    end

    def rdf_comment(v)
      statement(@primary_subject, RDF::RDFS.comment, RDF::Literal.new(v))
    end

    def rdf_ch_name(vs)
      statements = []
      vs.each do |v|
        if @ch_bnode == nil
          @ch_bnode = RDF::Node.new
          statements << statement(@primary_subject, @mbo[:compound], @ch_bnode)
          statements << statement(@ch_bnode, RDF.type, @mbo[:SampleChemicalCompound])
        end
        statements << statement(@ch_bnode, @mbo[:ch_name], RDF::Literal.new(v))
      end
      statements
    end

    def rdf_ch_compound_class(v)
      statement(@ch_bnode, @mbo[:ch_compound_class], RDF::Literal.new(v))
    end

    def rdf_ch_formula(v)
      statement(@ch_bnode, @mbo[:ch_formula], RDF::Literal.new(v))
    end

    def rdf_ch_exact_mass(v)
      statements = []
      bnode = RDF::Node.new
      statements << statement(@ch_bnode, @mbo[:ch_exact_mass], bnode)
      statements << statement(bnode, RDF.value,
                              RDF::Literal.new(v, :datatype => RDF::XSD.decimal))
      statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000088])
      statements
    end

    def rdf_ch_smiles(v)
      statement(@ch_bnode, @mbo[:ch_smiles], RDF::Literal.new(v))
    end

    def rdf_ch_iupac(v)
      statement(@ch_bnode, @mbo[:ch_iupac], RDF::Literal.new(v))
    end

    def rdf_ch_link(vs)
      statements = []
      vs.each do |v|
        /^(\S+)\s+(.+)$/ =~ v
        db = $1
        entry_id = $2
        objt = ""
        objt2 = ""
        case db
        when "CAS"
          objt = RDF::URI.new(@prefixes[:cas] + entry_id)
        when "CHEBI"
          objt = RDF::URI.new(@prefixes[:chebi] + entry_id)
        when "CHEMPDB"
          objt = RDF::URI.new(@prefixes[:chempdb] + entry_id)
        when "CHEMSPIDER"
          objt = RDF::URI.new(@prefixes[:chemspider] + entry_id)
        when "INCHIKEY"
          objt = RDF::URI.new(@prefixes[:inchikey] + entry_id)
        when "KEGG"
          objt = RDF::URI.new(@prefixes[:kegg] + entry_id)
        when "KNAPSACK"
          objt = RDF::URI.new(@prefixes[:knapsack] + entry_id)
        when "NIKKAJI"
          objt = RDF::URI.new(@prefixes[:nikkaji] + entry_id)
        when "LIPIDMAPS"
          objt = RDF::URI.new(@prefixes[:lipidmaps] + entry_id)
        when "LIPIDBANK"
          objt = RDF::URI.new(@prefixes[:lipidbank] + entry_id)
        when "PUBCHEM"
          if /CID:(\d+)\s+SID:(\d+)/ =~ entry_id
            objt = RDF::URI.new(@prefixes[:pubchem_compound] + $1)
            objt2 = RDF::URI.new(@prefixes[:pubchem] + $2)
          elsif /SID:(\d+)/ =~ entry_id
            objt = RDF::URI.new(@prefixes[:pubchem] + $1)
          elsif /CID:(\d+)/ =~ entry_id
            objt = RDF::URI.new(@prefixes[:pubchem_compound] + $1)
          else
            put_error_message("PUBCHEM CID or SID ? #{entry_id}")
          end
        when "HMDB"
          objt = RDF::URI.new(@prefixes[:hmdb] + entry_id)
        when "COMPTOX"
          objt = RDF::URI.new(@prefixes[:comptox] + entry_id)
        else
          put_error_message("Unlisted DB: #{db} .")
        end
        if objt2 != ""
          statements << statement(@ch_bnode, RDF::RDFS.seeAlso, objt)
          statements << statement(@ch_bnode, RDF::RDFS.seeAlso, objt2)
        elsif objt != ""
          statements << statement(@ch_bnode, RDF::RDFS.seeAlso, objt)
        else
          put_error_message("Could not found DB:#{db} ID:#{entry_id} .")
        end
      end
      statements
    end

    def rdf_sp_scientific_name(v)
      statements = []
      if @sp_bnode == nil
        @sp_bnode = RDF::Node.new
        statements << statement(@primary_subject, @mbo[:biological_sample], @sp_bnode)
        statements << statement(@sp_bnode, RDF.type, @mbo[:BiologicalSample])
        statements << statement(@sp_bnode, @mbo[:sp_scientific_name],
                                RDF::Literal.new(v, :language => :en))
      end
      statements
    end

    def rdf_sp_name(v)
      statements = []
      if @sp_bnode == nil
        @sp_bnode = RDF::Node.new
        statements << statement(@primary_subject, @mbo[:biological_sample], @sp_bnode)
        statements << statement(@sp_bnode, RDF.type, @mbo[:BiologicalSample])
        statements << statement(@sp_bnode, @mbo[:sp_name],
                                RDF::Literal.new(v, :language => :en))
        return statements
      else
        return statement(@sp_bnode, @mbo[:sp_name], RDF::Literal.new(v, :language => :en))
      end
    end

    def rdf_sp_lineage(v)
      statement(@sp_bnode, @mbo[:sp_lineage], RDF::Literal.new(v, :language => :en))
    end

    def rdf_sp_link(v)
      /(\S+)\s+(\d+)/ =~ v
      taxid = $2
      objt = RDF::URI.new(@prefixes[:taxonomy] + taxid)
      statement(@sp_bnode, RDF::RDFS.seeAlso, objt)
    end

    def rdf_sp_sample(v)
      statement(@sp_bnode, @mbo[:sp_sample], RDF::Literal.new(v, :language => :en))
    end

    def rdf_ac_instrument(v)
      statements = []
      if @ac_bnode == nil
        @ac_bnode = RDF::Node.new
        statements << statement(@primary_subject,
                                @mbo[:analytical_methods_and_conditions], @ac_bnode)
        statements << statement(@ac_bnode, RDF.type, @mbo[:AnalyticalMethods])
      end
      statements << statement(@ac_bnode, @mbo[:ac_instrument], RDF::Literal.new(v))
      statements
    end

    def rdf_ac_instrument_type(v)
      statement(@ac_bnode, @mbo[:instrument_type], RDF::Literal.new(v))
    end

    def rdf_ac_mass_spectrometry(vs)
      statements = []
      vs.each do |v|
        /^(\S+)\s+(.+)$/ =~ v
        subtag = $1
        value = $2
        case subtag
        when "MS_TYPE"
          statements << statement(@ac_bnode, @mbo[:ms_type], value)
        when "ION_MODE"
          statements << statement(@ac_bnode, @mbo[:ion_mode], value)
        when "COLLISION_ENERGY"
          unit = ""
          if /^(.+)\s+(\S+)$/ =~ value
            # ex. Ramp 21.1-31.6 eV
            collision_energy = $1
            unit = $2
            bnode = RDF::Node.new
            statements << statement(@ac_bnode, @mbo[:collision_energy], bnode)
            statements << statement(bnode, RDF.value,
                                    RDF::Literal.new(collision_energy))
          elsif /^(\d+)\s+(.+)$/ =~ value
            collision_energy = $1
            unit = $2
            bnode = RDF::Node.new
            statements << statement(@ac_bnode, @mbo[:collision_energy], bnode)
            statements << statement(bnode, RDF.value,
                                    RDF::Literal.new(collision_energy,
                                                     :datatype => RDF::XSD.decimal))
          end
          if unit == "kV"
            statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000248])
          elsif unit == "V"
            statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000218])
          elsif unit == "eV"
            statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000266])
          end

        when "COLLISION_GAS"
          statements << statement(@ac_bnode, @mbo[:collision_gas], value)
        when "DATE"
          statements << statement(@primary_subject, RDF::Vocab::DC.date,
                                  RDF::Literal.new(value.gsub(".", "-"),
                                                   :datatype => RDF::XSD.date))
        when "DESOLVATION_GAS_FLOW"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "DESOLVATION_TEMPERATURE"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "IONIZATION_ENERGY"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "LASER"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "MATRIX"
          statements << statement(@ac_bnode, @mbo[:matrix], value)
        when "MASS_ACCURACY"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "REAGENT_GAS"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "SCANNING"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "IONIZATION"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "IONIZATION_VOLTAGE"
          unit = ""
          if /^(.+)\s+(\S+)$/ =~ value
            # ex. 3.9 kV
            ionization_voltage = $1
            unit = $2
            bnode = RDF::Node.new
            statements << statement(@ac_bnode, @mbo[:ionization_voltage], bnode)
            statements << statement(bnode, RDF.value,
                                    RDF::Literal.new(ionization_voltage))
          elsif /^(\d+)\s+(.+)$/ =~ value
            ionization_voltage = $1
            unit = $2
            bnode = RDF::Node.new
            statements << statement(@ac_bnode, @mbo[:ionization_voltage], bnode)
            statements << statement(bnode, RDF.value,
                                    RDF::Literal.new(ionization_voltage,
                                                     :datatype => RDF::XSD.decimal))
          end
          if unit == "kV"
            statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000248])
          elsif unit == "V"
            statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000218])
          elsif unit == "eV"
            statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000266])
          end
        when "FRAGMENTATION_MODE"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        when "RESOLUTION"
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym],
                                  RDF::Literal.new(value, :datatype => RDF::XSD.integer))
        else
          put_error_message("Unknown subtag:#{subtag} .")
          statements << statement(@ac_bnode, @mbo[subtag.downcase.to_sym], value)
        end
      end
      statements
    end

    def rdf_ac_chromatography(vs)
      statements = []
      vs.each do |v|
        /^(\S+)\s+(.+)$/ =~ v
        subtag = $1
        dsc = $2
        case subtag
        when "CAPILLARY_VOLTAGE"
          /^([\d\.]+)\s+(\S+)$/ =~ dsc
          voltage = $1
          unit = $2
          bnode = RDF::Node.new
          statements << statement(@ac_bnode, @mbo[:capillary_voltage], bnode)
          statements << statement(bnode, RDF.value,
                                  RDF::Literal.new(voltage, :datatype => RDF::XSD.decimal))
          statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000248])
        when "COLUMN_NAME"
          statements << statement(@ac_bnode, @mbo[:column_name], RDF::Literal.new(dsc))
        when "COLUMN_TEMPERATURE"
          statements << statement(@ac_bnode, @mbo[:column_temperature], RDF::Literal.new(dsc))
        when "COLUMN_PRESSURE"
          statements << statement(@ac_bnode, @mbo[:column_pressure], RDF::Literal.new(dsc))
        when "FLOW_GRADIENT"
          statements << statement(@ac_bnode, @mbo[:flow_gradient], RDF::Literal.new(dsc))
        when "FLOW_RATE"
          /^([\d\.]+)\s+(\S+)$/ =~ dsc
          flow_rate = $1
          unit = $2
          bnode = RDF::Node.new
          statements << statement(@ac_bnode, @mbo[:flow_rate], bnode)
          statements << statement(bnode, RDF.value,
                                  RDF::Literal.new(flow_rate, :datatype => RDF::XSD.decimal))
          statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000271])
        when "RETENTION_TIME"
          /^([\d\.]+)\s+(\S+)$/ =~ dsc
          r_time = $1
          unit = $2
          bnode = RDF::Node.new
          statements << statement(@ac_bnode, @mbo[:retention_time], bnode)
          statements << statement(bnode, RDF.value,
                                  RDF::Literal.new(r_time, :datatype => RDF::XSD.decimal))
          statements << statement(bnode, @sio[:SIO_000221], @obo[:UO_0000031])
        when "SOLVENT"
          statements << statement(@ac_bnode, @mbo[:solvent], RDF::Literal.new(dsc))
        when "NAPS_RTI"
          statements << statement(@ac_bnode, @mbo[:naps_rti],
                                  RDF::Literal.new(dsc, :datatype => RDF::XSD.integer))
        else
          put_error_message("Unknown subtag:#{subtag} .")
        end
      end
      statements
    end

    def rdf_ms_focused_ion(vs)
      statements = []
      vs.each do |v|
        if @ms_bnode == nil
          @ms_bnode = RDF::Node.new
          statements << statement(@primary_subject, @mbo[:mass_spectram_meta_data], @ms_bnode)
          statements << statement(@ms_bnode, RDF.type, @mbo[:MassSpectramMetaData])
        end
        /^(\S+)\s+(.+)$/ =~ v
        subtag = $1
        value = $2
        case subtag
        when "BASE_PEAK"
          statements << statement(@ms_bnode, @mbo[subtag.downcase.to_sym],
                                  RDF::Literal.new(value, :datatype => RDF::XSD.decimal))
        when "DERIVATIVE_FORM"
          statements << statement(@ms_bnode, @mbo[subtag.downcase.to_sym], value)
        when "DERIVATIVE_MASS"
          statements << statement(@ms_bnode, @mbo[subtag.downcase.to_sym],
                                  RDF::Literal.new(value, :datatype => RDF::XSD.decimal))
        when "DERIVATIVE_TYPE"
          statements << statement(@ms_bnode, @mbo[subtag.downcase.to_sym], value)
        when "ION_TYPE"
          if @ion_type2id.key?(value)
            statements << statement(@ms_bnode, @mbo[:ion_type],
                                    @ion_type2id[value])
          else
            put_error_message("Unlisted ION_TYPE: #{value} .")
            statements << statement(@ms_bnode, @mbo[:ion_type_label],
                                    RDF::Literal.new(value))
          end
        when "PRECURSOR_M/Z"
          statements << statement(@ms_bnode, @mbo[:precursor_mz],
                                  RDF::Literal.new(value, :datatype => RDF::XSD.decimal))
        when "PRECURSOR_TYPE"
          if @ion_type2id.key?(value)
            statements << statement(@ms_bnode, @mbo[:precursor_type],
                                    @mbo[@ion_type2id[value].to_sym])
          else
            statements << statement(@ms_bnode, @mbo[:precursor_type_label],
                                    RDF::Literal.new(value))
            put_error_message("Unlisted PRECURSOR_TYPE: #{value} .")
          end
        else
          put_error_message("Unknown subtag:#{subtag} .")
        end
      end
      statements
    end

    def rdf_ms_data_processing(v)
      statement(@primary_subject, @mbo[:ms_data_processing],
                RDF::Literal.new(v, :language => :en))
    end

    def rdf_pk_splash(v)
      statement(@primary_subject, @mbo[:pk_splash], RDF::Literal.new(v))
    end

    def rdf_pk_annotation(vs)
      statements = []
      vs.each do |h|
        bnode = RDF::Node.new
        statements << statement(@primary_subject, @mbo[:has_peak_annotations], bnode)
        h.to_a.each do |e|
          case e[0]
          when /error\(ppm\)/
            bnode2 = RDF::Node.new
            statements << statement(bnode, @mbo[:error], bnode2)
            statements << statement(bnode2, @sio[:SIO_000221], @obo[:UO_0000169])
            statements << statement(bnode2, RDF.value,
                                    RDF::Literal(e[1], :datatype => RDF::XSD.decimal))
          when /mass/
            bnode2 = RDF::Node.new
            statements << statement(bnode, @mbo[:mass], bnode2)
            statements << statement(bnode2, @sio[:SIO_000221], @obo[:UO_0000221])
            statements << statement(bnode2, RDF.value,
                                    RDF::Literal(e[1], :datatype => RDF::XSD.decimal))
          when /m\/z/
            statements << statement(bnode, @mbo[:mz],
                                    RDF::Literal.new(e[1], :datatype => RDF::XSD.decimal))
          when /tentative_formula/
            statements << statement(bnode, @mbo[:tentative_formula], RDF::Literal.new(e[1]))
          when /formula_count/
            statements << statement(bnode, @mbo[:formula_count],
                                    RDF::Literal.new(e[1], :datatype => RDF::XSD.integer))
          when /num/
            statements << statement(bnode, @mbo[:num],
                                    RDF::Literal.new(e[1], :datatype => RDF::XSD.integer))
          else
            statements << statement(bnode, @mbo[e[0].to_sym], RDF::Literal.new(e[1]))
          end
        end
      end
      statements
    end

    def rdf_pk_num_peak(v)
      statement(@primary_subject, @mbo[:pk_num_peak],
                RDF::Literal.new(v, :datatype => RDF::XSD.integer))
    end

    def rdf_pk_peak(vs)
      statements = []
      vs.each_with_index do |v, i|
        /^(\S+)\s+(\S+)\s+(\S+)$/ =~ v
        mz = RDF::Literal.new($1, :datatype => RDF::XSD.decimal)
        int = RDF::Literal.new($2, :datatype => RDF::XSD.decimal)
        rel_int = RDF::Literal.new($3, :datatype => RDF::XSD.decimal)
        bnode = RDF::Node.new
        statements << statement(@primary_subject, @mbo[:has_peak], bnode)
        statements << statement(bnode, @mbo[:mz], mz)
        statements << statement(bnode, @mbo[:intensity], int)
        statements << statement(bnode, @mbo[:relative_intensity], rel_int)
        statements << statement(bnode, @mbo[:peak_sequential_number], RDF::Literal.new(i))
      end
      statements
    end

    def statement(subject, predicate, object)
      RDF::Statement.new(subject, predicate, object)
    end

    def put_error_message(msg)
      STDERR.puts "### [#{@accession}] :[#{@section}] [#{msg}]\n"
    end

  end

end
