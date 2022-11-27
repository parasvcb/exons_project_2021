from django.db import models


class Gene(models.Model):
    name = models.CharField(max_length=100, default="default")
    entrezid = models.IntegerField(default=0)
    organism = models.CharField(max_length=20, default="default")
    txid = models.IntegerField(default=0)


class Transcripts(models.Model):
    tId = models.CharField(max_length=20, default="default")
    swissprot = models.CharField(max_length=20, default="default")
    length = models.IntegerField(default=0)
    pi = models.BooleanField(default=0)
    exonscount = models.IntegerField(default=0)
    exonsIds = models.TextField(default="default")
    unique_domC = models.IntegerField(default=0)
    total_domC = models.IntegerField(default=0)
    exonscountUTMRD = models.CharField(max_length=100, default=0)
    exonscountAG = models.CharField(max_length=50, default=0)
    structured_count_disp = models.CharField(max_length=7, default=0)
    domains = models.TextField(default="default")
    exonsRegion = models.TextField(default="default")
    secondaryStructure = models.TextField(default="default")
    structured_count_ssp = models.CharField(max_length=7, default=0)
    geneRef = models.ForeignKey(
        Gene, related_name='trans', on_delete=models.CASCADE)


class ExonGenes(models.Model):
    wef = models.FloatField(default=0.0)
    exId = models.CharField(max_length=50, default="default")
    length = models.IntegerField(default=0)
    codst = models.IntegerField(default=0)
    codend = models.IntegerField(default=0)
    aaseq = models.TextField(default="default")
    rawst = models.IntegerField(default=0)
    rawend = models.IntegerField(default=0)
    sstype = models.IntegerField(default=0)
    parent = models.CharField(max_length=50, default="default")
    gene = models.ForeignKey(Gene, related_name='exonsGenes',
                             on_delete=models.CASCADE)


class ExonssPrediction(models.Model):
    list_trans_fk = models.TextField(default="default")
    ssseq = models.TextField(default="default")
    exon = models.ForeignKey(ExonGenes, related_name='exonsSS',
                             on_delete=models.CASCADE)


class ExonsDisorder(models.Model):
    list_trans_fk = models.TextField(default="default")
    disseq = models.TextField(default="default")
    exon = models.ForeignKey(ExonGenes, related_name='exonsDis',
                             on_delete=models.CASCADE)


class ExonsDomseq(models.Model):
    list_trans_fk = models.TextField(default="default")
    domseq = models.TextField(default="default")
    exon = models.ForeignKey(ExonGenes, related_name='exonsDom',
                             on_delete=models.CASCADE)


class DomainInfGene(models.Model):
    color = models.CharField(max_length=20, default="default")
    code = models.CharField(max_length=2, default="default")
    name = models.CharField(max_length=50, default="default")
    dId = models.CharField(max_length=50, default="default")
    gene = models.ForeignKey(
        Gene, related_name='domainCod', on_delete=models.CASCADE)
