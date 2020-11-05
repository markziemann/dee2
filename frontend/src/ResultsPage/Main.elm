module ResultsPage.Main exposing (..)

import Dict
import Maybe.Extra as MExtra
import ResultsPage.Helpers exposing (stageResultForDownload)
import ResultsPage.Types exposing (..)
import SearchPage.Helpers exposing (delay)
import Table


init : Model
init =
    { searchHits = Nothing
    , searchResultRows = Nothing
    , resultsTableQuery = ""
    , resultsTableState = Table.initialSort "id"
    , selectedResultsTableQuery = ""
    , selectedResultsTableState = Table.initialSort "id"
    , downloading = False
    , selectedResults = Dict.empty
    , paginationOffset =
        { perPage = 20
        , offset = 0
        }
    }


onlyData : Model -> ( Model, Cmd msg, Maybe OutMsg )
onlyData model =
    ( model, Cmd.none, Nothing )


update : Msg -> Model -> ( Model, Cmd Msg, Maybe OutMsg )
update msg model =
    case msg of
        ResultClicked result ->
            if Dict.member result.id model.selectedResults then
                onlyData
                    { model
                        | selectedResults = Dict.remove result.id model.selectedResults
                    }

            else
                onlyData
                    { model
                        | selectedResults =
                            MExtra.unwrap model.selectedResults
                                (\value -> Dict.insert result.id value model.selectedResults)
                                (stageResultForDownload result)
                    }

        SetResultsTableQuery resultsTableQuery ->
            onlyData { model | resultsTableQuery = resultsTableQuery }

        SetResultsTableState resultsTableState ->
            onlyData { model | resultsTableState = resultsTableState }

        SetSelectedResultsTableQuery selectedResultsTableQuery ->
            onlyData { model | selectedResultsTableQuery = selectedResultsTableQuery }

        SetSelectedResultsTableState selectedResultsTableState ->
            onlyData { model | selectedResultsTableState = selectedResultsTableState }

        DownloadRequested ->
            ( { model | downloading = True }, delay 5000 DownloadButtonReset, Nothing )

        DownloadButtonReset ->
            onlyData { model | downloading = False }
