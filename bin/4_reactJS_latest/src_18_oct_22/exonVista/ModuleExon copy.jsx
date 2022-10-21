import React, { Component, Fragment } from "react";
import ProtvistaManager from "protvista-manager";
import ProtvistaTrack from "protvista-track";
import ProtvistaNavigation from "protvista-navigation";
import ProtvistaSequence from "protvista-sequence";
import ProtvistaColouredSequence from "protvista-coloured-sequence";
import ProtvistaInterproTrack from "protvista-interpro-track";
import loadWebComponent from "../utils/load-web-component";
import sequence from "../mocks/sequence.json";
import { dataIPR, signatures, withResidues } from "../mocks/interpro";
import secondaryStructureData from "../mocks/interpro-secondary-structure.json";
import ProtvistaSaver from "protvista-saver";
import ProtvistaOverlay from "protvista-overlay";
import ProtvistaZoomTool from "protvista-zoom-tool";

// console.log('inmODULE')
class ModuleExon extends Component {
  componentDidMount() {
    document.querySelector("#interpro-track").data = dataIPR;
    document.querySelector("#interpro-track").contributors = signatures;
    document.querySelector("#interpro-track").fixedHighlight = "400:600";
    //
    document.querySelector("#sequence-track").data = sequence;
    document.querySelector("#sequence-track").fixedHighlight = "400:600";
    //
    document.querySelector("#sequence-coloured-track").data = sequence;
    document.querySelector("#sequence-coloured-track").fixedHighlight = "400:600";
    //
    document.querySelector("#sequence-coloured-track-iso").data = sequence;
    document.querySelector("#sequence-coloured-track-iso").fixedHighlight = "400:600";
    //
    document.querySelector("#track3").data = secondaryStructureData;
    document.querySelector("#track3").fixedHighlight = "400:600";

    //Includes a title in the exported file.
    document.querySelector("#saver").preSave = () => {
      const base = document.querySelector("#example");
      const title = document.createElement("h2");
      title.setAttribute("id", "tmp_title_element");
      title.innerHTML = "ProtVista Snapshot";
      console.log("Paras")
      console.log("-->",title,base.firstChild)
      base.insertBefore(title, base.firstChild);
    };
    //removes the title from the DOM
    document.querySelector("#saver").postSave = () => {
      document
        .querySelector("#example")
        .removeChild(document.getElementById("tmp_title_element"));
    };

    //Sets the background color of the image to save.
    document.querySelector("#saver").backgroundColor = "#ffffff";
    document.querySelector("#saver2").backgroundColor = "#ddddee";
  }
  //feeding data into properties, 
  //
  
  render() {
    // console.log('out100',document.querySelector("#track3"))
    loadWebComponent("protvista-manager", ProtvistaManager);
    loadWebComponent("protvista-track", ProtvistaTrack);
    loadWebComponent("protvista-navigation", ProtvistaNavigation);
    loadWebComponent("protvista-sequence", ProtvistaSequence);
    loadWebComponent("protvista-coloured-sequence", ProtvistaColouredSequence);
    loadWebComponent("protvista-interpro-track", ProtvistaInterproTrack);
    loadWebComponent("protvista-saver", ProtvistaSaver);
    loadWebComponent("protvista-overlay", ProtvistaOverlay);
    loadWebComponent("protvista-zoom-tool", ProtvistaZoomTool);
    //console.log('out200',document.querySelector("#interpro-track"))
    return (
      <Fragment>
        <protvista-saver element-id="example" id="saver" />
        <protvista-saver
          element-id="just-tracks"
          id="saver2"
          file-name="tracks"
          file-format="jpeg"
        >
          {/* {console.log('out300')} */}
          <button>Download Just Tracks</button>
        </protvista-saver>
        <protvista-overlay for="just-tracks" />
        <protvista-manager
        //length = ?? can be added
          //attributes="variantfilters"
          attributes="paras"
          //displaystart="-10"
          displaystart="10"
          displayend="100"
          id="example"
        >
          {/* {console.log('out400')} */}
          <protvista-zoom-tool
            length="790"
            style={{
              float: "right",
              "--button-background": "#00639a",
              "--button-text-color": "#FFFFFF",
              "--button-background-focus": "#00a6d5",
              "--button-border-radius": "4px",
            }}
          >
            {/* {console.log('out500')} */}
            <span slot="zoom-in">+</span>
            <span slot="zoom-in-seq">Zoom to Sequence</span>
          </protvista-zoom-tool>
          <protvista-navigation length="790" />
          {/* length="456"
          displaystart="143"
          displayend="400"
          highlightStart="23"
          highlightEnd="45"
          rulerstart="50", can be added */}
          <div id="just-tracks">
            <protvista-sequence
              length="790"
              id="sequence-track"
              highlight-event="onmouseover"
              use-ctrl-to-zoom
            />
            <protvista-coloured-sequence
              length="790"
              id="sequence-coloured-track"
              scale="hydrophobicity-interface-scale"
              height="10"
              highlight-event="onmouseover"
              use-ctrl-to-zoom
            />
            <protvista-coloured-sequence
              length="790"
              id="sequence-coloured-track-iso"
              scale="isoelectric-point-scale"
              color_range="white:0,dodgerblue:11"
              height="10"
              use-ctrl-to-zoom
            />
            <protvista-track
            //https://github.com/ebi-webcomponents/nightingale/tree/master/packages/protvista-track
            // have nice components, shape, data to be loaded
              id="track3"
              length="790"
              displaystart="1"
              displayend="80"
              layout="non-overlapping"
              use-ctrl-to-zoom
            />
            <protvista-interpro-track
            // start: number (optional)
            // The start position of the selected region.
            // end: number (optional)
              id="interpro-track"
              length="790"
              shape="roundRectangle"
              highlight-event="onmouseover"
              expanded
              use-ctrl-to-zoom
            />
          </div>
        </protvista-manager>
      </Fragment>
    );
  }
}

export default ModuleExon;
