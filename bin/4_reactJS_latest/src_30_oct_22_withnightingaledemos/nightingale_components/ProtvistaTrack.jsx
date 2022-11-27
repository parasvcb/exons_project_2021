import React, { Fragment, Component } from "react";
import DataLoader from "dataloader";
import ProtvistaFeatureAdapter from "protvista-feature-adapter";
import ProtvistaTrack from "protvista-track";
import loadWebComponent from "../utils/load-web-component";
import data from "../mocks/features.json";
import data2 from "../mocks/features2.json";
import data3 from "../mocks/interpro-secondary-structure.json";

window.customElements.define("data-loader", DataLoader);
//window.customElements.define("protvista-track", ProtvistaTrack);

console.log("logD",document)
class ProtvistaTrackWrapper extends Component {
  componentDidMount() {
    console.log('0')
    document.querySelector("#track1").data = data;
    document.querySelector("#track2").data = data;
    document.querySelector("#track3").data = data;
    document.querySelector("#track4").data = data2;
    document.querySelector("#track5").data = data3;
  }

  render() {
    console.log('p2')
    loadWebComponent("protvista-feature-adapter", ProtvistaFeatureAdapter);
    loadWebComponent("protvista-track", ProtvistaTrack);
    loadWebComponent("data-loader", DataLoader);
    console.log("p1")
    return (
      <Fragment>
        <h2>Track with data-loader</h2>
        <protvista-track
          length="770"
          displaystart="1"
          displayend="770"
          tooltip-event="click"
        >
          <protvista-feature-adapter id="adapter1">
            <data-loader>
              <source src="https://www.ebi.ac.uk/proteins/api/features/P05067?categories=PTM" />
            </data-loader>
          </protvista-feature-adapter>
        </protvista-track>
        <h2>Track with set data</h2>
        <h3>Default layout</h3>
        <protvista-track id="track1" width="1020" length="223" />
        <h3>Non-overlapping layout</h3>
        <protvista-track id="track2" length="223" layout="non-overlapping" />
        <h3>Multiple highlights</h3>
        <protvista-track id="track3" length="223" highlight="20:30,40:80" />
        <h3>Shapes</h3>
        <protvista-track
          id="track4"
          length="100"
          layout="non-overlapping"
          shape="roundRectangle"
          height="420"
        />
        <h3>Track with InterPro data model</h3>
        <protvista-track
          id="track5"
          length="100"
          width="1020"
          height="100"
          displaystart="1"
          displayend="80"
          layout="non-overlapping"
        />
      </Fragment>
    );
  }
}

export default ProtvistaTrackWrapper;
