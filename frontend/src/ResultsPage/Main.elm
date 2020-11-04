module ResultsPage.Main exposing (..)

import Array
import ResultsPage.Types exposing (..)
import SearchPage.Helpers exposing (delay)
import Table


init : Model
init =
    { searchHits = Nothing
    , searchResultRows = Nothing
    , resultsTableState = Table.initialSort "id"
    , resultsTableQuery = ""
    , downloading = False
    }


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    case msg of
        ResultClicked result ->
            let
                searchResultRows =
                    Maybe.map (Array.set result.id { result | selected = not result.selected })
                        model.searchResultRows
            in
            ( { model | searchResultRows = searchResultRows }, Cmd.none )

        SetResultsTableQuery resultsTableQuery ->
            ( { model | resultsTableQuery = resultsTableQuery }
            , Cmd.none
            )

        SetResultsTableState resultsTableState ->
            ( { model | resultsTableState = resultsTableState }
            , Cmd.none
            )

        DownloadRequested ->
            ( { model | downloading = True }, delay 5000 DownloadButtonReset )

        DownloadButtonReset ->
            ( { model | downloading = False }, Cmd.none )
