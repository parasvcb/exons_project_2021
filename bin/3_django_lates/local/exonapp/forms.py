from django import forms



class GeneInpForm(forms.Form):
    NCBI_gene_id = forms.IntegerField(required=False)
    Gene_name = forms.CharField(max_length=100, required=False)

    def clean(self):
        #cleaned_data = super().clean()
        cleaned_data=super(GeneInpForm,self).clean()
        gennam = cleaned_data.get("Gene_name")
        entid = cleaned_data.get("NCBI_gene_id")

        if entid and gennam:
            # only one should be entered
            raise forms.ValidationError("Please enter either entrez id or genename not both")

