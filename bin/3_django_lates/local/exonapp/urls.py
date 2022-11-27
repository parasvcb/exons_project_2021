import datetime
from . import views
from django.urls import path
print("\n\nURLS\n\n")
now = datetime.datetime.now()

#handler404 = 'exonapp.views.page_not_found'
#handler500 = 'exonapp.views.server_error'
#handler403 = 'exonapp.views.permission_denied'

urlpatterns = [
    path('gene/<int:pk>/', views.DetailGene.as_view(), name='Gene-detail'),
    path('name/<str:name>/', views.GeneName.as_view(), name='Gene-namesearch'),
    path('ncbid/<int:entrezid>/', views.GeneIdRet.as_view(),
         name='Gene-ncbidsearch'),
    path('', views.catchall)
]
print("urls Passed")
print(now)
    #path('stats/', views.Stats.as_view())
