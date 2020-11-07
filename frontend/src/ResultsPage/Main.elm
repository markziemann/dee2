module ResultsPage.Main exposing (..)

import Dict
import Dict.Extra as DExtra
import Maybe.Extra as MExtra
import ResultsPage.Helpers exposing (stageResultForDownload)
import ResultsPage.Types exposing (..)
import SearchPage.Helpers exposing (delay)
import Set
import Table
import SharedTypes exposing (RemoteData(..))

init : Model
init =
    { searchResults = NotAsked
    , resultsTableQuery = ""
    , resultsTableState = Table.initialSort "id"
    , selectedResultsTableQuery = ""
    , selectedResultsTableState = Table.initialSort "id"
    , downloading = False
    , selectedResults = Dict.empty
    , resultsPendingRemoval = Set.empty
    , paginationOffset =
        { perPage = 20
        , offset = 0
        }
    }


onlyData : Model -> ( Model, Cmd msg, OutMsg)
onlyData model =
    ( model, Cmd.none, Nothing)


update : Msg -> Model -> ( Model, Cmd Msg, OutMsg)
update msg model =
    case msg of
        ResultClicked result ->
            if Dict.member result.id model.selectedResults then
                onlyData
                    { model
                        | selectedResults = Dict.remove result.id model.selectedResults
                        , resultsPendingRemoval = Set.remove result.id model.resultsPendingRemoval
                    }

            else
                onlyData
                    { model
                        | selectedResults =
                            MExtra.unwrap model.selectedResults
                                (\value -> Dict.insert result.id value model.selectedResults)
                                (stageResultForDownload result)
                    }

        SelectedResultClicked id ->
            if Set.member id model.resultsPendingRemoval then
                onlyData { model | resultsPendingRemoval = Set.remove id model.resultsPendingRemoval }

            else
                onlyData { model | resultsPendingRemoval = Set.insert id model.resultsPendingRemoval }

        RemoveStagedSelections ->
            onlyData
                { model
                    | selectedResults =
                        DExtra.removeMany
                            (Set.intersect
                                (Set.fromList <| Dict.keys model.selectedResults)
                                model.resultsPendingRemoval
                            )
                            model.selectedResults
                    , resultsPendingRemoval = Set.empty
                }

        PageRequest paginationOffset ->
            ( model, Cmd.none, Just paginationOffset)

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
