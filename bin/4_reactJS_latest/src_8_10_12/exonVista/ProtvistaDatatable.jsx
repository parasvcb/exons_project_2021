import React, { useEffect, useState } from "react";
import ProtvistaDatatable from "protvista-datatable";
import ProtvistaTrack from "protvista-track";
import ProtvistaManager from "protvista-manager";
import ProtvistaNavigation from "protvista-navigation";
import { transformData } from "protvista-feature-adapter";
import { load } from "data-loader";
import loadWebComponent from "../utils/load-web-component";
import { Fragment } from "react";
import { data } from "../mocks/datatable";
const ProtvistaDatatableWrapper = () => {

  loadWebComponent("protvista-datatable", ProtvistaDatatable);
  loadWebComponent("protvista-manager", ProtvistaManager);
  loadWebComponent("protvista-track", ProtvistaTrack);
  loadWebComponent("protvista-navigation", ProtvistaNavigation);
   console.log("Inside")
  // useEffect(() => {
  //   async function fetchData() {
  //     const loadedData = await load(
  //       "https://www.ebi.ac.uk/proteins/api/features/P05067?categories=MOLECULE_PROCESSING"
  //     );
  //     console.log("fresh",loadedData)
  //     const transformedData = transformData(loadedData?.payload);
  //     //remove first array element
  //     //transformedData=transformData.shift()
  //     console.log("data", transformedData)
  //     setData(transformedData);
  //   }
  //   fetchData();
  // }, []);



  return (
    <>
      <h2>Track with data-loader</h2>
      <protvista-manager attributes="length displaystart displayend highlight selectedid">
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
              {console.log("dataNew",data)}
              {data?.map((row, i) => (
                <Fragment key={i}>
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
                    <td>Something hidden</td>
                  </tr>
                </Fragment>
              ))}
            </tbody>
          </table>
        </protvista-datatable>
      </protvista-manager>
    </>
  );
};

export default ProtvistaDatatableWrapper;
