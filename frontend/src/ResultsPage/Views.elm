module ResultsPage.Views exposing (..)

import Array
import Dict
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events as Events exposing (onClick)
import ResultsPage.Helpers exposing (queryString)
import ResultsPage.Types exposing (..)
import SearchPage.Types exposing (SearchResult, SearchResults)
import Table exposing (HtmlDetails, Status)


type alias SelectedResult =
    ( Int, ( String, String ) )


selectClickedResult : SelectedResults -> SearchResult -> List (Html.Attribute Msg)
selectClickedResult selectedResults result =
    [ class
        (if Dict.member result.id selectedResults then
            "table-primary"

         else
            ""
        )
    , style "cursor" "pointer" -- First place not using bootstrap for style?
    , Events.onClick (ResultClicked result)
    ]


get : String -> SearchResult -> String
get key searchResult =
    Maybe.withDefault "" (Dict.get key searchResult.data)


getId =
    .id >> String.fromInt


defaultTable =
    Table.defaultCustomizations


resultsTable : SelectedResults -> Table.Config SearchResult Msg
resultsTable selectedResults =
    Table.customConfig
        { toId = getId
        , toMsg = SetResultsTableState
        , columns =
            [ Table.stringColumn "Row" getId
            , Table.stringColumn "Species" (get "species")
            , Table.stringColumn "SRA Run" (get "SRR_accession")
            , Table.stringColumn "QC summary" (get "QC_summary")
            , Table.stringColumn "SRA experiment" (get "SRX_accession")
            , Table.stringColumn "SRA sample" (get "SRS_accession")
            , Table.stringColumn "SRA project" (get "SRP_accession")
            , Table.stringColumn "Sample" (get "Sample_name")
            , Table.stringColumn "Experiment" (get "GEO_series")

            --GEO_series
            ]
        , customizations =
            { defaultTable
                | tableAttrs =
                    [ class "table table-hover table-sm table-bordered table-responsive-xl"
                    , style "table-layout" "fixed" -- Prevents table going wider than parent element
                    ]
                , rowAttrs = selectClickedResult selectedResults
            }
        }


selectedTable : Table.Config SelectedResult Msg
selectedTable =
    Table.customConfig
        { toId = Tuple.first >> String.fromInt
        , toMsg = SetSelectedResultsTableState
        , columns =
            [ Table.intColumn "Row" Tuple.first
            , Table.stringColumn "Species" (Tuple.second >> Tuple.first)
            , Table.stringColumn "SRA" (Tuple.second >> Tuple.second)
            ]
        , customizations =
            { defaultTable
                | tableAttrs =
                    [ class "table table-sm table-bordered"
                    ]
            }
        }


hideWhenTrue : String -> Bool -> String
hideWhenTrue classString true =
    if true then
        String.join " " [ classString, "invisible" ]

    else
        classString


buttonOrSpinner : Bool -> SelectedResults -> Html Msg
buttonOrSpinner downloading rows =
    if not downloading then
        a
            [ hideWhenTrue "btn btn-outline-primary btn-block" (Dict.isEmpty rows) |> class
            , href <| "/download/" ++ queryString rows
            , attribute "download" "data.zip"
            , onClick DownloadRequested
            ]
            [ text "download" ]

    else
        button
            [ class "btn btn-primary btn-block"
            , attribute "type" "button"
            , attribute "disabled" "disabled"
            ]
            [ span
                [ class "spinner-border spinner-border-sm mr-2"
                , attribute "role" "status"
                , attribute "aria-hidden" "true"
                ]
                []
            , text "Loading..."
            ]


viewSearchResults : Model -> List (Html Msg)
viewSearchResults ({ searchResultRows, resultsTableState, resultsTableQuery, searchHits } as model) =
    [ div [ class "row" ]
        [ div [ class "col-xl-10" ]
            [ div [ class "bg-light text-primary" ]
                [ text "Hits: ", text (Maybe.withDefault "" (Maybe.map String.fromInt searchHits)) ]
            , Table.view
                (resultsTable model.selectedResults)
                resultsTableState
                (Maybe.withDefault [] (Maybe.map Array.toList searchResultRows))
            ]
        , div [ class "col-xl-2" ]
            [ div [class "sticky-top"]
                [ div [ class "bg-light text-primary" ]
                    [ text <| "Selected: " ++ (Dict.size model.selectedResults |> String.fromInt) ]
                , Table.view
                    selectedTable
                    model.selectedResultsTableState
                    (Dict.toList model.selectedResults)
                , div [ class "btn-group mb-4" ] [ buttonOrSpinner model.downloading model.selectedResults ]
                ]
            ]
        ]
    ]
