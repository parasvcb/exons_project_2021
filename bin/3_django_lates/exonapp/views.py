from cgi import print_environ
from django.http import HttpResponse, HttpResponseRedirect, Http404, HttpResponseForbidden
from django.views.generic import TemplateView
from django.template import RequestContext
from django.shortcuts import render
from .forms import *
from .models import *
from rest_framework import generics
from .serializers import *
import datetime
print("\n\nViews\n\n")


#from django.core.urlresolvers import reverse
#from django.shortcuts import render_to_response, redirect, render


# Create your views here.


class DetailGene(generics.RetrieveAPIView):
    queryset = Gene.objects.all()
    serializer_class = GeneSerializer


catchall = TemplateView.as_view(template_name='index.html')


class CatchAll(TemplateView):
    template_name = "index.html"


class GeneIdRet(generics.RetrieveAPIView):
    print('**\n\n**')
    # print(self)
    queryset = Gene.objects.all()
    serializer_class = GeneSerializer
    # print_environ
    lookup_field = 'entrezid'
# when above is called, then queryset is built, means it shuld call all the fields of Gene.objects.all(), it needs output informat of geneSerializer and match with lookup field

class GeneName(generics.ListAPIView):
    queryset = Gene.objects.all()
    serializer_class = GeneSerializerNames
    now = datetime.datetime.now()
    print ("\n\nOverHere\n\n")
    print (now)
    def get_queryset(self):
        print("selfkwargs", self.kwargs)
        gpat = self.kwargs['name']
        print ("gpat", gpat)
        # return Gene.objects.all()
        print(Gene.objects.all())
        return self.queryset.filter(name__icontains=gpat)[:200]
# we should query the queryset which is all of the gene objects, then it should present results in form of the geneSerializernames, but what it should query isnt clear, 
# it wueries in def queryset(), 

# class Stats(generics.ListAPIView):
#     queryset = Gene.objects.all()
#     serializer_class = GeneSerializerNames

#     def get_queryset(self):

#         print("selfkwargs", self.kwargs)
#         gpat = self.kwargs['pk']
#         # return Gene.objects.all()
#         print(Gene.objects.all())
#         return self.queryset.filter(name__icontains=gpat)[:200]
## can be worked on 