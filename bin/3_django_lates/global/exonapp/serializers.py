from rest_framework import serializers
from .models import *

print("\nserializers\n")


class DomainInfGeneSerializer(serializers.ModelSerializer):
    class Meta:
        fields = ('color', 'code', 'name', 'dId')
        model = DomainInfGene


class TranscriptSerializer(serializers.ModelSerializer):
    class Meta:
        # add fieldfs of transcipt object only
        fields = ('tId', 'id', 'swissprot', 'length', 'pi', 'structured_count_disp', 'structured_count_ssp',
                  'exonsIds', 'exonscount', 'unique_domC', 'total_domC', 'exonscountUTMRD', 'exonscountAG', 'domains','exonsRegion','secondaryStructure')
        model = Transcripts


class ExonssPredictionSerializer(serializers.ModelSerializer):
    class Meta:
        fields = ('list_trans_fk', 'ssseq')
        model = ExonssPrediction


class ExonsDisorderSerializer(serializers.ModelSerializer):
    class Meta:
        fields = ('list_trans_fk', 'disseq')
        model = ExonsDisorder


class ExonsDomseqSerializer(serializers.ModelSerializer):
    class Meta:
        fields = ('list_trans_fk', 'domseq')
        model = ExonsDomseq


class ExonGeneSerializer(serializers.ModelSerializer):
    exonsDom = ExonsDomseqSerializer(many=True)
    exonsDis = ExonsDisorderSerializer(many=True)
    exonsSS = ExonssPredictionSerializer(many=True)

    class Meta:
        fields = ('exId', 'length', 'parent', 'aaseq', 'wef', 'id', 'exonsDom',
                  'exonsDis', 'exonsSS', 'codst', 'codend', 'rawst',
                  'rawend', 'sstype')

        #fields = ('sstype')

        model = ExonGenes


class GeneSerializer(serializers.ModelSerializer):
    trans = TranscriptSerializer(many=True)
    exonsGenes = ExonGeneSerializer(many=True)
    domainCod = DomainInfGeneSerializer(many=True)

    class Meta:
        fields = (
            'name', 'id', 'organism', 'entrezid', 'txid', 'trans', 'exonsGenes', 'domainCod'
        )
        model = Gene


class GeneSerializerNames(serializers.ModelSerializer):
    class Meta:
        fields = (
            'name', 'id', 'organism', 'txid', 'entrezid'
        )
        model = Gene


