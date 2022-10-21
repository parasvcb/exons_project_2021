import React, { Component } from 'react'
import ProtvistaTrackWrapper from '../nightingale_components/ProtvistaTrack'

export class Track extends Component {
  render() {
    console.log("in")
    return (
        <div>
        <p>?IsIt?</p>
        <ProtvistaTrackWrapper/>
        <p>?IsIt?</p>
        </div>
    )
    }
}

export default Track