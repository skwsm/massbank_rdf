#!/usr/bin/env ruby
require 'rubygems'
require 'rdf'
require 'rdf/ntriples'
require 'rdf/turtle'
require '../lib/massbank.rb'
include RDF

f = open(ARGV.shift)
f_out = open(ARGV.shift, "w")
r = MassBank::Record.new(f.read)
#r.record.each do |e|
#  p e
#end
factory = MassBank::RDFFactory.new(r.record, f_out)
factory.rdfize

