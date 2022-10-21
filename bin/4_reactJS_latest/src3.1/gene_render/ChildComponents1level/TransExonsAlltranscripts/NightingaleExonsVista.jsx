import React, { Component, Fragment, useState, useEffect, useCallback } from "react";
import ProtvistaManager from "protvista-manager";
import ProtvistaTrack from "protvista-track";
import ProtvistaNavigation from "protvista-navigation";
import ProtvistaSequence from "protvista-sequence";
import ProtvistaInterproTrack from "protvista-interpro-track";
import ProtvistaColouredSequence from "protvista-coloured-sequence";
import loadWebComponent from "./utils/load-web-component";
//import sequence from "../mocks/sequence.json";
// import { dataIPR, signatures } from "../mocks/interpro";
// import secondaryStructureData from "../mocks/interpro-secondary-structure.json";
import ProtvistaSaver from "protvista-saver";
import ProtvistaOverlay from "protvista-overlay";
import ProtvistaZoomTool from "protvista-zoom-tool";
import ProtvistaDatatable from "protvista-datatable";
//import { dataTable } from "../mocks/datatable";

/*
Three tracks, exons, domains, secondary structres and disorder, exonTrack
*/
class NightingaleExonsVista extends Component {
  constructor(props){
    super(props);
    //console.log ('props',this.props)
    this.state = {
      sequence: this.props.sequence,
      domains: this.props.nightDom,
      ss: this.props.nightSs,
      exons: this.props.nighthTable,
  }}
  componentDidMount() {
     console.log("cdm1",this.props.transwhole)
     console.log("cdm2",this.props.transwhole.nightDom)
     //this["interpro-track"].innerHTML = this.props.transwhole.nightDom;
     document.querySelector("#sequence-track").data = this.props.transwhole.sequence;
     document.querySelector("#domain-track").data = this.props.transwhole.nightDom;
     //document.querySelector("#exon-track").data = this.props.transwhole.nightExons;
     document.querySelector("#secondaryStructure-track").data = this.props.transwhole.nightSs;
     //document.querySelector("#datatable").data = this.props.transwhole.nightExons
     
  }
    //document.querySelector("#interpro-track").data = this.props.nightDom;
    //* document.querySelector("#interpro-track").contributors = signatures;
    //* document.querySelector("#interpro-track").fixedHighlight = "460:480, 500:520";
    //
    //* document.querySelector("#sequence-track").data = sequence;
    //document.querySelector("#sequence-track").fixedHighlight = "400:600";
    //
    //document.querySelector("#sequence-coloured-track").data = sequence;
    //document.querySelector("#sequence-coloured-track").fixedHighlight = "400:600";
    //
    //document.querySelector("#sequence-coloured-track-iso").data = sequence;
    //document.querySelector("#sequence-coloured-track-iso").fixedHighlight = "400:600";
    //
    //* document.querySelector("#track3").data = secondaryStructureData;
    //document.querySelector("#track3").fixedHighlight = "400:600";

    //Includes a title in the exported file.


    // document.querySelector("#saver").preSave = () => {
    //   console.log ('*No')
    //   const base = document.querySelector("#example");
    //   const title = document.createElement("h2");
    //   title.setAttribute("id", "tmp_title_element");
    //   title.innerHTML = "ProtVista Snapshot";
    //   console.log("Paras")
    //   console.log("-->",title,base.firstChild)
    //   base.insertBefore(title, base.firstChild);
    // };
    
    
    //removes the title from the DOM


    // document.querySelector("#saver").postSave = () => {
    //   document
    //     .querySelector("#example")
    //     .removeChild(document.getElementById("tmp_title_element"));
    // };


    //Sets the background color of the image to save.

    // document.querySelector("#saver").backgroundColor = "#ffffff";
    // document.querySelector("#saver2").backgroundColor = "#ddddee";
  // }
  render () {
    loadWebComponent("protvista-manager", ProtvistaManager);
    //loadWebComponent("protvista-table", ProtvistaDatatableWrapper);
    loadWebComponent("protvista-track", ProtvistaTrack);
    loadWebComponent("protvista-navigation", ProtvistaNavigation);
    loadWebComponent("protvista-sequence", ProtvistaSequence);
    //loadWebComponent("protvista-coloured-sequence", ProtvistaColouredSequence);
    loadWebComponent("protvista-interpro-track", ProtvistaInterproTrack);
    loadWebComponent("protvista-saver", ProtvistaSaver);
    loadWebComponent("protvista-overlay", ProtvistaOverlay);
    loadWebComponent("protvista-zoom-tool", ProtvistaZoomTool);
    loadWebComponent("protvista-datatable", ProtvistaDatatable);
    //document.querySelector("#interpro-track").data = this.props.nightDom
    console.log("render",this.props.nightDom)
    console.log("render2",this.props)
    const lengthSeq=this.props.transwhole.sequence.length
    return (
      <Fragment>
         <protvista-manager
        //length = ?? can be added
          //attributes="variantfilters"
          //attributes="paras" 
          displaystart= "0" displayend="100" length ={lengthSeq} //id="example"
          //displaystart="-10"
          // displaystart="10"
          // displayend="100"
          // id="example"
        >
          <protvista-navigation length ={lengthSeq}  /> 

          <div id="just-tracks">
          <protvista-sequence
            length={lengthSeq} displaystart="0" displayend="100"
            id="sequence-track"
            highlight-event="onmouseover"
            use-ctrl-to-zoom
          />
           <protvista-track
              id="secondaryStructure-track"
              length={lengthSeq}
              layout="non-overlapping"
              use-ctrl-to-zoom
            />
            <protvista-interpro-track
              id="domain-track"
              length={lengthSeq}
              highlight-event="onmouseover"
              expanded
              use-ctrl-to-zoom
            />
          </div>
        </protvista-manager>

      </Fragment>
    )
  }
  }

export default NightingaleExonsVista;
