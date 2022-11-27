import React, { Component, Fragment, useState, useEffect, useCallback } from "react";
import ProtvistaManager from "protvista-manager";
import ProtvistaTrack from "protvista-track";
import ProtvistaNavigation from "protvista-navigation";
import ProtvistaSequence from "protvista-sequence";
import ProtvistaColouredSequence from "protvista-coloured-sequence";
import ProtvistaInterproTrack from "protvista-interpro-track";
import loadWebComponent from "../utils/load-web-component";
import sequence from "../mocks/sequence.json";
import { dataIPR, signatures } from "../mocks/interpro";
import secondaryStructureData from "../mocks/interpro-secondary-structure.json";
import ProtvistaSaver from "protvista-saver";
import ProtvistaOverlay from "protvista-overlay";
import ProtvistaZoomTool from "protvista-zoom-tool";
import ProtvistaDatatable from "protvista-datatable";
import { dataTable } from "../mocks/datatable";

//import ProtvistaDatatableWrapper from "./ProtvistaDatatable";

class ModuleExon extends Component {
  componentDidMount() {
    //document.querySelector("#linegraph").data = linegraphData;

    document.querySelector("#interpro-track").data = dataIPR;
    document.querySelector("#interpro-track").contributors = signatures;
    document.querySelector("#interpro-track").fixedHighlight = "460:480, 500:520";
    //
    document.querySelector("#sequence-track").data = sequence;
    //document.querySelector("#sequence-track").fixedHighlight = "400:600";
    //
    //document.querySelector("#sequence-coloured-track").data = sequence;
    //document.querySelector("#sequence-coloured-track").fixedHighlight = "400:600";
    //
    //document.querySelector("#sequence-coloured-track-iso").data = sequence;
    //document.querySelector("#sequence-coloured-track-iso").fixedHighlight = "400:600";
    //
    document.querySelector("#track3").data = secondaryStructureData;
    //document.querySelector("#track3").fixedHighlight = "400:600";

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
  render() {
    // console.log('out100',document.querySelector("#track3"))
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
    //loadWebComponent("nightingale-linegraph-track", NightingaleLinegraphTrack);
    //loadWebComponent("nightingale-linegraph-track", NightingaleLinegraphTrack);
    

    //console.log('out200',document.querySelector("#interpro-track"))
    return (
      <Fragment>
        <protvista-saver element-id="example" id="saver" />
        <protvista-saver
          element-id="just-tracks" id="saver2" file-name="tracks" file-format="jpeg"
        >
          <button>Download Just Tracks</button>
        </protvista-saver>
        {/* <protvista-overlay for="just-tracks" /> */}
        <protvista-manager
        //length = ?? can be added
          //attributes="variantfilters"
          //attributes="paras" 
          displaystart= "1" displayend="400" length ="870" //id="example"
          //displaystart="-10"
          // displaystart="10"
          // displayend="100"
          // id="example"
        >
          <protvista-zoom-tool
            length="870"
            displaystart="50"
            displayend="770"
            style={{
              float: "right",
              "--button-background": "#00639a",
              "--button-text-color": "#FFFFFF",
              "--button-background-focus": "#00a6d5",
              "--button-border-radius": "4px",
            }}
          >
            <span slot="zoom-in">+</span>
            <span slot="zoom-in-seq">Zoom to Sequence</span>
          </protvista-zoom-tool>
          <protvista-navigation length="870"  /> 

          <div id="just-tracks">
            <protvista-sequence
              length="870" displaystart="0" displayend="770"
              id="sequence-track"
              highlight-event="onmouseover"
              use-ctrl-to-zoom
            />
            {/* <protvista-coloured-sequence
              length="870"
              id="sequence-coloured-track"
              scale="hydrophobicity-interface-scale"
              height="10"
              start="-50"
              end="90"
              highlight-event="onmouseover"
              use-ctrl-to-zoom
            />
            <protvista-coloured-sequence
              length="870"
              start="-50"
              end="90"
              id="sequence-coloured-track-iso"
              scale="isoelectric-point-scale"
              color_range="white:0,dodgerblue:11"
              height="10"
              use-ctrl-to-zoom
            /> */}
            <protvista-track
              id="track3"
              length="870"
              layout="non-overlapping"
              use-ctrl-to-zoom
            />
            <protvista-interpro-track
              id="interpro-track"
              length="870"
              shape="roundRectangle"
              layout="non-overlapping"
              highlight-event="onmouseover"
              expanded
              use-ctrl-to-zoom
            />
          </div>
          {/* <ProtvistaDatatableWrapper/> */}
          {/* <protvista-table/> */}
        
        <protvista-datatable displaystart="-50" displayend="820">
          <table>
            <thead>
              <tr>
                <th data-filter="ft_key">Feature key</th>
                <th>Description</th>
                <th>Positions</th>
              </tr>
            </thead>
            <tbody>
              {/* {console.log("dataNew",dataTable)} */}
              {dataTable?.map((row, i) => (
                //i is row index, starting from 0
                <Fragment key={i}>
                  {console.log("da*",row,i)}
                  <tr
                    data-id={`${row.start}-${row.end}`}
                    data-start={row.start}
                    data-end={row.end}
                  >
                    <td data-filter="ft_key" data-filter-value={row.type}>
                      {row.type}
                    </td>
                    <td>{row.description}</td>
                    <td>
                      {row.start}-{row.end}
                    </td>
                  </tr>
                  <tr data-group-for={`${row.start}-${row.end}`}>
                    <td>{row.tooltipContent}</td>
                  </tr>
                </Fragment>
              ))}
            </tbody>
          </table>
        </protvista-datatable>
        </protvista-manager>
      </Fragment>
      
    );
  }
}

export default ModuleExon;

          /* navigation attributes 
          length="456"
          displaystart="143"
          displayend="400"
          highlightStart="23"
          highlightEnd="45"
          rulerstart="50", can be added 

          track in protvista
          //https://github.com/ebi-webcomponents/nightingale/tree/master/packages/protvista-track
          // have nice components, shape, data to be loaded
          displaystart="-50"
          displayend="90"
              
          interpro track
            // start: number (optional)
            // The start position of the selected region.
            // end: number (optional)
            //tooltip="onmouseover"
              tooltip-event="click"
        
          <nightingale-linegraph-track id="linegraph" length="870" />
        */