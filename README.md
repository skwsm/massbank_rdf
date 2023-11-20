# massbank_rdf

Example usage

```
> git clone https://github.com/skwsm/massbank_rdf.git

> git clone https://github.com/MassBank/MassBank-data.git
```

chdir to the lib directory

```
> for i in ../../MassBank-data/*/*.txt; do name=`basename -s txt $i`; ruby ../bin/massbank2rdf.rb $i ../output/${name}ttl; done
```

